package vcfcmd

import (
	"bufio"
	"encoding/json"
	"fmt"
	"os"
	"sort"
	"strings"

	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/compgenlab/cgkit/internal/varstore"
	"github.com/spf13/cobra"
)

var (
	vcfVarQueryOutput   string
	vcfVarQuerySamples  []string
	vcfVarQueryVariants []string
	vcfVarQueryRegion   string
	vcfVarQueryMinDP    int
	vcfVarQueryHomRef   bool
	vcfVarQueryFormat   string
	vcfVarQueryStore    string
	vcfVarQueryVerbose  bool
)

var vcfVarQueryCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": "v0.5.0"},
	Use:         "vcf-varquery <input.vcf | store-base>",
	Short:       "Query which subjects carry a variant, or which variants a subject carries",
	Long: `Query genotypes without caring which format holds them. The input may be
a VCF (plain or bgzipped) or a Parquet store written by vcf-toparquet. A store
may be named by its base ("cohort" for cohort.calls.parquet...), by its
directory ("cohort/" or "cohort" for cohort/calls.parquet...), or by any one of
its three member files. The backend is inferred from the path; override with
--store.

Two modes, one of which must be given:

  --sample NAME     report the variants that subject carries
  --variant LOCUS   report the subjects carrying that variant (chrom:pos:ref:alt)

Both are repeatable. --min-dp gates calls by depth; a gate over data lacking DP
admits everything rather than rejecting it, since an absent depth is not evidence
of a shallow one.

There is deliberately no --min-gq. GQ is recorded per ALT call, but conversion
builds its callable runs from depth alone, so no GQ survives for a reference call
-- a store could not honor the gate there while a VCF would, and the two backends
would silently disagree. The gq column is in the output, so filter on it
downstream if you need to.

By default only the alternate-allele calls are reported -- "which variants does
this subject carry". --hom-ref switches to every interrogated site -- "show me
all the sites for this subject" -- adding the reference calls to the same
stream: in --variant mode the samples that are 0/0 at the locus, in --sample mode
the sites that subject was 0/0 at. The gt column tells the two apart.

A 0/0 row means the whole genotype was reference. At a multiallelic record a 0/2
sample is not a carrier of allele 1, yet it is not reference either, so it
appears under neither -- writing 0/0 for it would be a genotype the source never
contained.

In --sample mode this walks the entire sites catalog, since being reference is a
statement about every site the callset interrogated -- expect roughly one row per
variant in the source. Use --region to bound it.

A reference call is only as good as the evidence for it: the gate must admit it
(a 0/0 at DP 3 under --min-dp 10 is not a reference observation), an off-catalog
locus yields nothing at all, and an incomplete Parquet store refuses rather than
guesses -- reporting an unobserved sample as reference would invent an
observation. One asymmetry is unavoidable: a store keeps only alternate
genotypes, so a reference call recovered from one carries a synthesized 0/0 with
no DP/AD/GQ, where the same query against a VCF reports the recorded genotype
and its quality fields.

A VCF can always answer, because it stores an explicit genotype for every sample
at every record. A Parquet store reconstructs the reference calls from its sites
and callable-regions files, and refuses if either is missing or it was built
without coverage information.

  --sample NAME       subject to report variants for (repeatable)
  --variant LOCUS     variant to report carriers of (repeatable)
  --region R          restrict to a 1-based region (chrom:start-end)
  --min-dp N          minimum depth for a call to count
  --hom-ref           report every interrogated site, not only the ALT calls
  --format F          tsv (default), json, vcf, or list
  --store KIND        force the backend: vcf or parquet
  -v, --verbose       report the backend, the store's conversion settings, and
                      whether the quality gate could actually act, on stderr

Both query modes emit one layout, so --sample and --variant output can be
concatenated, sorted and cut the same way:

  chrom pos ref alt sample gt dp min_dp ad_ref ad_alt gq

The locus leads as four columns rather than one packed chrom:pos:ref:alt field,
so rows cut and sort on position directly. The chromosome is echoed the way it
was asked for, whichever naming convention the file itself uses. By default only
valid ALT calls are reported: genotypes carrying the alternate and passing the
gate.

min_dp is the tightest lower bound on depth the backend can vouch for. Where a
call records its own depth that is the bound. A reference call recovered from a
store has no recorded depth, but it came from a run built at the conversion
--min-dp, so the bound is that threshold -- which is the evidence that the site
was covered well enough to call. Anything else reports ".", since a threshold
written where nothing is known would assert a depth the data never had.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		nSample, nVariant := len(vcfVarQuerySamples), len(vcfVarQueryVariants)
		if nSample == 0 && nVariant == 0 {
			return fmt.Errorf("give at least one --sample or --variant")
		}
		if nSample > 0 && nVariant > 0 {
			return fmt.Errorf("--sample and --variant are separate modes; use one at a time")
		}
		format := strings.ToLower(vcfVarQueryFormat)
		switch format {
		case "tsv", "json":
		case "vcf":
			if nVariant > 0 {
				return fmt.Errorf("--format vcf applies to --sample mode, which emits variants")
			}
		case "list":
			if nSample > 0 {
				return fmt.Errorf("--format list applies to --variant mode, which emits sample ids")
			}
			// A bare list of ids cannot say which of them are carriers and which
			// are reference calls, and a mixed list read as carriers is exactly the
			// conflation the two categories exist to prevent.
			if vcfVarQueryHomRef {
				return fmt.Errorf("--format list emits bare sample ids, which cannot distinguish " +
					"a carrier from a reference call; use --format tsv or json with --hom-ref")
			}
		default:
			return fmt.Errorf("unknown format %q (use tsv, json, vcf, or list)", format)
		}

		store, err := openVarStore(args[0], vcfVarQueryStore)
		if err != nil {
			return err
		}
		defer store.Close()

		span, err := varstore.ParseSpan(vcfVarQueryRegion)
		if err != nil {
			return err
		}
		gate := varstore.Gate{MinDP: int32(vcfVarQueryMinDP)}
		if vcfVarQueryVerbose {
			describeStore(cmd, store, args[0], gate)
		}

		w, closeFn, err := openOutput(cmd, vcfVarQueryOutput)
		if err != nil {
			return err
		}
		out := bufio.NewWriter(w)

		if nSample > 0 {
			err = runVarQuerySamples(out, store, span, gate, format, args[0])
		} else {
			warnUnknownSites(cmd, store)
			err = runVarQueryVariants(cmd, out, store, gate, format, args[0])
		}
		if err != nil {
			return err
		}

		if err := out.Flush(); err != nil {
			return err
		}
		if closeFn != nil {
			return closeFn()
		}
		return nil
	},
}

// openVarStore picks a backend for path. kind may force "vcf" or "parquet";
// empty infers from the filename.
func openVarStore(path, kind string) (varstore.Store, error) {
	switch strings.ToLower(kind) {
	case "vcf":
		return varstore.OpenVcf(path)
	case "parquet":
		return varstore.OpenParquet(path)
	case "":
	default:
		return nil, fmt.Errorf("unknown store %q (use vcf or parquet)", kind)
	}

	// A directory-form store: "cohort/", or the directory itself.
	if varstore.IsDirBase(path) {
		return varstore.OpenParquet(path)
	}
	if st, err := os.Stat(path); err == nil && st.IsDir() {
		return varstore.OpenParquet(path)
	}

	lower := strings.ToLower(path)
	switch {
	case strings.HasSuffix(lower, ".parquet"):
		return varstore.OpenParquet(path)
	case strings.HasSuffix(lower, ".vcf"), strings.HasSuffix(lower, ".vcf.gz"),
		strings.HasSuffix(lower, ".vcf.bgz"), strings.HasSuffix(lower, ".bcf"):
		return varstore.OpenVcf(path)
	}
	// A bare base name is a Parquet store when its calls file exists.
	if _, err := os.Stat(varstore.CallsPath(varstore.TrimStoreSuffix(path))); err == nil {
		return varstore.OpenParquet(path)
	}
	return nil, fmt.Errorf("cannot tell what kind of store %q is; pass --store vcf or --store parquet", path)
}

// provenance writes the ## header lines shared by the tabular outputs.
func provenance(out *bufio.Writer, source string) {
	fmt.Fprintln(out, "## program: "+buildinfo.String())
	fmt.Fprintln(out, "## cmd: "+buildinfo.CommandLine())
	fmt.Fprintln(out, "## input: "+source)
}

// describeStore reports which backend was chosen and, for a Parquet store, the
// conversion settings baked into it.
//
// The --min-dp comparison is the point of this. A store's called-site runs were
// built at one threshold; a query gating at a different one is not asking a
// question the store can answer consistently, and the two backends will stop
// agreeing. That is invisible without being told.
func describeStore(cmd *cobra.Command, store varstore.Store, path string, g varstore.Gate) {
	out := cmd.ErrOrStderr()
	switch s := store.(type) {
	case *varstore.ParquetStore:
		p := s.Provenance()
		fmt.Fprintf(out, "store    parquet %s (%d samples)\n", path, p.NumSamples)
		if p.Source != "" {
			fmt.Fprintf(out, "  built from  %s\n", p.Source)
		}
		if p.Program != "" {
			fmt.Fprintf(out, "  by          %s\n", p.Program)
		}
		fmt.Fprintf(out, "  spans       %s", p.Spans)
		if p.Spans == varstore.SpansSites {
			fmt.Fprintf(out, "  (answers only for variants in the sites catalog)")
		}
		fmt.Fprintln(out)
		if p.NoCallable {
			fmt.Fprintf(out, "  callable    not tracked (--no-callable): --hom-ref will refuse\n")
		} else {
			fmt.Fprintf(out, "  min-dp      %d at conversion\n", p.MinDP)
			if g.MinDP > 0 && g.MinDP != p.MinDP {
				fmt.Fprintf(out, "  NOTE: querying at --min-dp %d but the runs were built at %d;\n"+
					"        non-carrier vs not-assayed will not match a direct VCF query.\n",
					g.MinDP, p.MinDP)
			}
		}
		// The one place the backends cannot agree, so it is worth saying rather
		// than leaving a column of dots to be puzzled over.
		if vcfVarQueryHomRef {
			fmt.Fprintf(out, "  NOTE: --hom-ref rows report a synthesized 0/0 with no DP/AD/GQ;\n"+
				"        a store keeps only alternate genotypes, so the reference call\n"+
				"        itself was never recorded -- only the fact that it was made.\n")
		}
	case *varstore.VcfStore:
		names, _ := s.Samples()
		fmt.Fprintf(out, "store    vcf %s (%d samples)\n", path, len(names))
	}
}

// reportGate says whether the quality gate could actually act.
//
// A gate over a field the data lacks admits everything rather than rejecting it,
// so a --min-dp that looks like a filter can be doing nothing at all. This is the
// one diagnostic most worth surfacing: the numbers look plausible either way, and
// only the field census distinguishes them.
func reportGate(cmd *cobra.Command, g varstore.Gate, calls []varstore.Call, excluded int) {
	if g.IsZero() {
		return
	}
	out := cmd.ErrOrStderr()
	var haveDP int
	for _, c := range calls {
		if c.DP != varstore.Missing {
			haveDP++
		}
	}
	fmt.Fprintf(out, "gate     excluded %d call(s)\n", excluded)
	if g.MinDP > 0 && haveDP == 0 && len(calls) > 0 {
		fmt.Fprintf(out, "  WARNING: --min-dp %d had no effect; no call here carries DP\n", g.MinDP)
	}
}

// warnUnknownSites notes any queried variant the source never reported.
//
// Without this, such a variant returns zero carriers -- or all not-assayed
// under --hom-ref -- which reads exactly like a real negative result. The
// source simply never looked there, and only a gVCF could say otherwise.
func warnUnknownSites(cmd *cobra.Command, store varstore.Store) {
	for _, v := range vcfVarQueryVariants {
		locus, err := varstore.ParseLocus(v)
		if err != nil {
			continue // reported properly by the query itself
		}
		known, err := store.SiteKnown(locus)
		if err != nil || known {
			continue
		}
		fmt.Fprintf(cmd.ErrOrStderr(),
			"warning: %s is not in the source; reporting not-assayed for every sample.\n"+
				"         A VCF only supports queries for the variants it contains.\n", locus)
	}
}

// varQuerySampleResult is one subject's carried variants.
type varQuerySampleResult struct {
	Sample string          `json:"sample"`
	Calls  []varstore.Call `json:"calls"`
}

// runVarQuerySamples reports the variants carried by each requested subject.
func runVarQuerySamples(out *bufio.Writer, store varstore.Store, span *varstore.Span,
	gate varstore.Gate, format, source string) error {

	var results []varQuerySampleResult

	for _, s := range vcfVarQuerySamples {
		q := varstore.Query{
			Samples:    []string{s},
			Gate:       gate,
			IncludeRef: vcfVarQueryHomRef,
		}
		if span != nil {
			q.Spans = []varstore.Span{*span}
		}
		// No sort: the store emits in its own (chrom, pos, alt, sample) order,
		// which is contig order rather than lexicographic.
		calls, err := varstore.CollectCalls(store, q)
		if err != nil {
			return err
		}
		results = append(results, varQuerySampleResult{Sample: s, Calls: calls})
	}

	switch format {
	case "json":
		enc := json.NewEncoder(out)
		enc.SetIndent("", "  ")
		return enc.Encode(results)
	case "vcf":
		return writeVarQueryVcf(out, results, source)
	}

	provenance(out, source)
	fmt.Fprintln(out, strings.Join(tsvColumns, "\t"))
	minDP := vouchedMinDP(store)
	for _, r := range results {
		for _, c := range r.Calls {
			fmt.Fprintf(out, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				c.Chrom, c.Pos, c.Ref, c.Alt, c.SampleID, c.GT,
				intOrDot(c.DP), minDPOf(c, minDP),
				intOrDot(c.ADRef), intOrDot(c.ADAlt), intOrDot(c.GQ))
		}
	}
	return nil
}

// writeVarQueryVcf emits a minimal VCF of the calls reported for each sample.
//
// It is deliberately spare. Without --hom-ref only carried variants are present,
// because a reference genotype was not asked for and must not be invented: a
// locus the requested samples do not carry is simply absent, never written as
// 0/0. With --hom-ref the reference calls are in the result set, and appear as
// the genotypes they are. A sample with no call at a locus another sample brought
// into the output stays "./.", which is the honest reading either way.
func writeVarQueryVcf(out *bufio.Writer, results []varQuerySampleResult, source string) error {
	names := make([]string, 0, len(results))
	for _, r := range results {
		names = append(names, r.Sample)
	}
	fmt.Fprintln(out, "##fileformat=VCFv4.2")
	fmt.Fprintln(out, "##source="+buildinfo.String())
	fmt.Fprintln(out, "##cgkit_vcf-varqueryCommand="+buildinfo.CommandLine())
	fmt.Fprintln(out, "##cgkit_vcf-varquerySource="+source)
	fmt.Fprintln(out, `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`)
	fmt.Fprintln(out, `##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">`)
	fmt.Fprintln(out, `##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">`)
	fmt.Fprintln(out, "#"+strings.Join(append([]string{
		"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
	}, names...), "\t"))

	// index calls by locus so one row can carry every requested sample
	type key struct {
		chrom string
		pos   int32
		ref   string
		alt   string
	}
	order := []key{}
	byLocus := map[key]map[string]varstore.Call{}
	for _, r := range results {
		for _, c := range r.Calls {
			k := key{c.Chrom, c.Pos, c.Ref, c.Alt}
			if _, ok := byLocus[k]; !ok {
				byLocus[k] = map[string]varstore.Call{}
				order = append(order, k)
			}
			byLocus[k][c.SampleID] = c
		}
	}
	// Sort rather than emit in first-seen order. Loci are gathered per sample, so
	// the second sample's private loci would otherwise all land after the first
	// sample's -- an unsorted VCF, which cannot be indexed.
	//
	// NOTE: this compares chromosomes LEXICOGRAPHICALLY, so chr10 sorts before
	// chr2. Harmless while a store holds few contigs, wrong for a real one, since
	// an indexable VCF needs the header's contig order. The fix is to emit in the
	// store's own order -- which Calls already streams in -- rather than to patch
	// this comparator, which cannot know the contig order.
	sort.Slice(order, func(i, j int) bool {
		a, b := order[i], order[j]
		if a.chrom != b.chrom {
			return a.chrom < b.chrom
		}
		if a.pos != b.pos {
			return a.pos < b.pos
		}
		if a.ref != b.ref {
			return a.ref < b.ref
		}
		return a.alt < b.alt
	})
	for _, k := range order {
		cols := []string{k.chrom, fmt.Sprint(k.pos), ".", k.ref, k.alt, ".", ".", ".", "GT:DP:GQ"}
		for _, n := range names {
			if c, ok := byLocus[k][n]; ok {
				cols = append(cols, fmt.Sprintf("%s:%s:%s", c.GT, intOrDot(c.DP), intOrDot(c.GQ)))
			} else {
				cols = append(cols, "./.:.:.")
			}
		}
		fmt.Fprintln(out, strings.Join(cols, "\t"))
	}
	return nil
}

// runVarQueryVariants reports which subjects carry each requested variant.
func runVarQueryVariants(cmd *cobra.Command, out *bufio.Writer, store varstore.Store,
	gate varstore.Gate, format, source string) error {

	type variantResult struct {
		Variant string          `json:"variant"`
		Locus   varstore.Locus  `json:"locus"`
		Calls   []varstore.Call `json:"calls,omitempty"`
	}
	var results []variantResult

	for _, v := range vcfVarQueryVariants {
		locus, err := varstore.ParseLocus(v)
		if err != nil {
			return err
		}
		r := variantResult{Variant: locus.String(), Locus: locus}
		at := varstore.Query{Loci: []varstore.Locus{locus}, Gate: gate}

		// Compare against the ungated result so the gate's effect is a measured
		// number rather than an inference from the row count.
		var ungated []varstore.Call
		if vcfVarQueryVerbose && !gate.IsZero() {
			ungated, err = varstore.CollectCalls(store, varstore.Query{Loci: at.Loci})
			if err != nil {
				return err
			}
		}

		at.IncludeRef = vcfVarQueryHomRef
		calls, err := varstore.CollectCalls(store, at)
		if err != nil {
			return err
		}
		r.Calls = calls

		// Counted apart from the row total: the gate's effect is measured against
		// the ALT calls alone, and reference rows are not gated calls that
		// survived, they are a different question's answer.
		nCarriers, nHomRef := 0, -1
		if vcfVarQueryHomRef {
			nHomRef = 0
		}
		for _, c := range calls {
			if varstore.IsAltCarrier(c.GT) {
				nCarriers++
			} else if nHomRef >= 0 {
				nHomRef++
			}
		}
		if vcfVarQueryVerbose {
			reportVariant(cmd, store, locus, nCarriers, nHomRef)
			if !gate.IsZero() {
				reportGate(cmd, gate, ungated, len(ungated)-nCarriers)
			}
		}
		results = append(results, r)
	}

	switch format {
	case "json":
		enc := json.NewEncoder(out)
		enc.SetIndent("", "  ")
		return enc.Encode(results)
	case "list":
		seen := map[string]bool{}
		for _, r := range results {
			for _, c := range r.Calls {
				if !seen[c.SampleID] {
					seen[c.SampleID] = true
					fmt.Fprintln(out, c.SampleID)
				}
			}
		}
		return nil
	}

	provenance(out, source)
	minDP := vouchedMinDP(store)
	fmt.Fprintln(out, strings.Join(tsvColumns, "\t"))
	for _, r := range results {
		for _, c := range r.Calls {
			// The queried locus, not the call's stored spelling: every row of a
			// --variant query then reads back the way it was asked for, and both
			// sub-modes agree even though only this one has calls to draw on.
			fmt.Fprintf(out, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				r.Locus.Chrom, r.Locus.Pos, r.Locus.Ref, r.Locus.Alt,
				c.SampleID, c.GT,
				intOrDot(c.DP), minDPOf(c, minDP),
				intOrDot(c.ADRef), intOrDot(c.ADAlt), intOrDot(c.GQ))
		}
	}
	return nil
}

// intOrDot renders a possibly-missing integer field.
func intOrDot(v int32) string {
	if v == varstore.Missing {
		return "."
	}
	return fmt.Sprint(v)
}

// tsvColumns is the one tabular layout the genotype outputs use. The locus
// leads, as four columns rather than one packed chrom:pos:ref:alt field, so rows
// sort and cut on position without re-splitting a composite key.
var tsvColumns = []string{
	"chrom", "pos", "ref", "alt", "sample", "gt", "dp", "min_dp", "ad_ref", "ad_alt", "gq",
}

// vouchedMinDP is the depth threshold a store's callable runs were built at, or 0
// for a backend with no such threshold -- a VCF, where a call either carries its
// own depth or none is knowable.
func vouchedMinDP(store varstore.Store) int32 {
	if ps, ok := store.(*varstore.ParquetStore); ok {
		if p := ps.Provenance(); !p.NoCallable {
			return p.MinDP
		}
	}
	return 0
}

// minDPOf renders min_dp: the tightest lower bound on depth the backend can
// vouch for at this call.
//
// An exact depth is its own best bound. Where there is none, a store can still
// vouch for the threshold its runs were built at -- but only for a reconstructed
// reference call, which by construction came from such a run. A store's calls
// file holds only ALT-carrying genotypes, so an all-reference GT identifies those
// rows unambiguously. Anything else without a depth reports ".", because nothing
// is known: putting a threshold there would assert a depth the data never had.
func minDPOf(c varstore.Call, storeMinDP int32) string {
	if c.DP != varstore.Missing {
		return fmt.Sprint(c.DP)
	}
	if storeMinDP > 0 && c.GT == varstore.HomRefGT {
		return fmt.Sprint(storeMinDP)
	}
	return "."
}

func init() {
	f := vcfVarQueryCmd.Flags()
	f.StringVarP(&vcfVarQueryOutput, "output", "o", "-", "Output filename (- for stdout)")
	f.StringArrayVar(&vcfVarQuerySamples, "sample", nil, "Report variants carried by this subject (repeatable)")
	f.StringArrayVar(&vcfVarQueryVariants, "variant", nil, "Report subjects carrying this variant, as chrom:pos:ref:alt (repeatable)")
	f.StringVar(&vcfVarQueryRegion, "region", "", "Only variants in this 1-based region (chrom:start-end, or chrom); requires a tabix-indexed file")
	f.IntVar(&vcfVarQueryMinDP, "min-dp", 0, "Minimum DP for a call to count")
	f.BoolVar(&vcfVarQueryHomRef, "hom-ref", false, "Also report reference (0/0) calls, not only alternate carriers")
	f.StringVar(&vcfVarQueryFormat, "format", "tsv", "Output format: tsv, json, vcf, or list")
	f.StringVar(&vcfVarQueryStore, "store", "", "Force the backend: vcf or parquet (default: infer from the path)")
	f.BoolVarP(&vcfVarQueryVerbose, "verbose", "v", false, "Report the backend, store provenance and gate effect on stderr")
}

// reportVariant summarises one queried locus: what the catalog knows about it,
// how many ALT calls were reported, and how many reference rows --hom-ref
// contributed (homRef is -1 when it was not asked for).
func reportVariant(cmd *cobra.Command, store varstore.Store, l varstore.Locus,
	nCarriers, homRef int) {

	out := cmd.ErrOrStderr()
	fmt.Fprintf(out, "variant  %s\n", l)

	// Allele frequency comes from the catalog, which knows the site even when
	// nobody carries it -- so this stays meaningful at AC 0.
	if ps, ok := store.(*varstore.ParquetStore); ok {
		if site, found, err := ps.Site(l); err == nil && found {
			fmt.Fprintf(out, "  site        AC=%d AN=%d AF=%.6g  n_carriers=%d n_called=%d n_lowdp=%d\n",
				site.AC, site.AN, site.AF(), site.NCarriers, site.NCalled, site.NLowDP)
		}
	}
	fmt.Fprintf(out, "  alt calls   %d\n", nCarriers)
	if homRef >= 0 {
		fmt.Fprintf(out, "  hom-ref     %d\n", homRef)
	}
}
