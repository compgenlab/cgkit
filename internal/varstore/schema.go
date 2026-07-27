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
//	BASE.regions.parquet   runs of catalog sites at which a sample was called
//
// # What a store built from a plain VCF can answer
//
// Only the variants the VCF actually contains. A plain VCF asserts nothing
// about the positions between its records: an absent base was not observed to
// be reference, it simply was not reported. So the sites catalog is the exact
// boundary of what is knowable, and a query for any locus outside it yields
// StateNotAssayed for every sample -- not StateNonCarrier.
//
// This is why the run intervals in the regions file must never be read as
// coverage. They compress a per-sample, per-site record of "this sample was
// called at these catalog sites"; the gaps between those sites are not part of
// the claim. A gVCF is different, because its reference blocks (END, MIN_DP)
// are positive statements about spans, and only such a store could answer
// off-catalog positions.
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

// CalledSiteRun records that one sample was successfully called, at adequate
// depth, at every catalog site in [Start, End]. Start and End are the first and
// last such site positions, and NSites is how many catalog sites the run covers.
//
// IT SAYS NOTHING ABOUT THE GAPS BETWEEN THOSE SITES. A plain VCF records
// nothing between its variant records, so there is no basis for calling an
// intervening base reference -- the caller may simply never have looked. The
// interval form is a COMPRESSION of a per-sample, per-site fact ("called here,
// and here, and here"), not a statement about genomic territory. Reading it as
// a coverage interval would silently manufacture reference observations for
// positions that were never interrogated, which is the precise error the
// four-way classification exists to prevent.
//
// Consequently a run is only meaningful at positions that appear in the sites
// catalog, and Classify checks the catalog first for exactly that reason. A
// gVCF, whose reference blocks carry END and MIN_DP, is what would license
// answering off-catalog positions; see SpanSemantics.
type CalledSiteRun struct {
	SampleID string `parquet:"sample_id,dict"`
	Chrom    string `parquet:"chrom,dict"`
	Start    int32  `parquet:"start"`
	End      int32  `parquet:"end"`
	NSites   int32  `parquet:"n_sites"`
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
