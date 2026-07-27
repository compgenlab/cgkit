// Package varstore holds the on-disk schema for a sparse genotype store and a
// uniform way to query genotypes regardless of which format backs them.
//
// A Parquet store is a set of three files derived from one base name. All three
// are required: the calls file alone cannot distinguish a confidently-called
// reference genotype from a position that was never assayed, which is the
// distinction Classify exists to make.
//
//	BASE.calls.parquet     one row per ALT-carrying genotype
//	BASE.sites.parquet     one row per interrogated site, sample-independent
//	BASE.regions.parquet   one row per contiguous run of adequately-covered sites
//
// The sites file cannot be reconstructed from the calls file. Taking the
// distinct loci out of the calls recovers every site only when the store holds
// the entire joint callset; over a subset of samples the sites where nobody in
// that subset carries an ALT vanish silently, and a query would then report
// "never interrogated" for a position that was in fact interrogated.
package varstore

// Missing marks an absent integer field (DP, GQ, AD) in a Parquet row. VCFs
// routinely omit these -- a GT-only phased panel has no DP at all -- and the
// columns are kept non-optional so reads stay a flat scan, so the absence has
// to be in-band. Callers must test against Missing before using a value; a
// naive comparison would treat it as an extremely low quality score.
const Missing int32 = -1

// Call is one ALT-carrying genotype for one sample at one biallelic site.
// Records are split so that every row carries exactly one ALT allele.
type Call struct {
	SampleID string `parquet:"sample_id,dict"`
	Chrom    string `parquet:"chrom,dict"`
	Pos      int32  `parquet:"pos"`
	Ref      string `parquet:"ref,dict"`
	Alt      string `parquet:"alt,dict"`
	GT       string `parquet:"gt,dict"`
	DP       int32  `parquet:"dp"`
	ADRef    int32  `parquet:"ad_ref"`
	ADAlt    int32  `parquet:"ad_alt"`
	GQ       int32  `parquet:"gq"`
}

// Locus returns the site this call belongs to.
func (c Call) Locus() Locus {
	return Locus{Chrom: c.Chrom, Pos: c.Pos, Ref: c.Ref, Alt: c.Alt}
}

// Site is one interrogated variant site, independent of any sample. NCarriers,
// NLowDP and NCalled are counted across every sample present in the source, so
// a site with NCarriers == 0 still records that the position was examined.
type Site struct {
	Chrom     string `parquet:"chrom,dict"`
	Pos       int32  `parquet:"pos"`
	Ref       string `parquet:"ref,dict"`
	Alt       string `parquet:"alt,dict"`
	NCarriers int32  `parquet:"n_carriers"`
	NLowDP    int32  `parquet:"n_lowdp"`
	NCalled   int32  `parquet:"n_called"`
}

// Locus returns this site's identity.
func (s Site) Locus() Locus {
	return Locus{Chrom: s.Chrom, Pos: s.Pos, Ref: s.Ref, Alt: s.Alt}
}

// CallableRegion is a maximal run of consecutive variant sites at which one
// sample had DP at or above the conversion threshold. Start and End are the
// first and last such site positions.
//
// The span between two in-run sites is assumed callable, not observed: a
// plain VCF records nothing between its variant records. Only a gVCF, with
// reference blocks carrying END and MIN_DP, would make these intervals
// observations rather than interpolations.
type CallableRegion struct {
	SampleID string `parquet:"sample_id,dict"`
	Chrom    string `parquet:"chrom,dict"`
	Start    int32  `parquet:"start"`
	End      int32  `parquet:"end"`
	NSites   int32  `parquet:"n_sites"`
}

// Covers reports whether pos falls inside this run.
func (r CallableRegion) Covers(chrom string, pos int32) bool {
	return r.Chrom == chrom && pos >= r.Start && pos <= r.End
}

// File suffixes of the three members of a Parquet store.
const (
	CallsSuffix   = ".calls.parquet"
	SitesSuffix   = ".sites.parquet"
	RegionsSuffix = ".regions.parquet"
)

// CallsPath returns the calls file for a store base name.
func CallsPath(base string) string { return base + CallsSuffix }

// SitesPath returns the sites file for a store base name.
func SitesPath(base string) string { return base + SitesSuffix }

// RegionsPath returns the callable-regions file for a store base name.
func RegionsPath(base string) string { return base + RegionsSuffix }
