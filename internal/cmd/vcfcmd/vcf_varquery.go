package vcfcmd

import (
	"bufio"
	"encoding/json"
	"fmt"
	"os"
	"strings"

	"github.com/compgenlab/cghts/varstore"
	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/spf13/cobra"
)

var (
	vcfVarQueryOutput   string
	vcfVarQuerySamples  []string
	vcfVarQueryVariants []string
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

Two independent axes, at least one of which must be given. --variant selects
sites and --sample selects subjects; they compose rather than exclude.

--variant takes any of these, repeatably, and there is no separate --region:

  chr1                     a whole contig
  chr1:1000-2000           a region
  chr1:1000                any variant at that position
  chr1:1000:A:T            one exact variant
  panel.vcf                a file -- its format is detected from the content

A value is a file when one exists by that name, and an inline selector otherwise,
so a mistyped locus still gets a locus error rather than "no such file". Three
file formats are recognised:

  VCF/BCF        announced by its ##fileformat line; one target per ALT allele
  BED            chrom start end, 0-BASED half-open, as BED always is
  site list      whitespace-separated "chrom pos [ref alt]", 1-based; a line
                 holding a single token is an inline selector, so a file may just
                 list the tokens above

BED and site lists are told apart by their third column: a BED end coordinate is
numeric, a REF allele never is. Both may carry extra columns, '#' comments and
blank lines. Since a misread would shift coordinates by one rather than fail, -v
reports which format each file was read as and how many targets it produced, and
a file with no targets is an error rather than an empty result.

--sample alone reports what those subjects carry; --variant alone reports who
carries those variants; both together ask what those subjects' genotypes are at
those sites. An axis left unnamed is unrestricted, not empty.

--min-dp gates calls by depth; a gate over data lacking DP admits everything
rather than rejecting it, since an absent depth is not evidence of a shallow one.

There is deliberately no --min-gq. GQ is recorded per ALT call, but conversion
builds its callable runs from depth alone, so no GQ survives for a reference call
-- a store could not honor the gate there while a VCF would, and the two backends
would silently disagree. The gq column is in the output, so filter on it
downstream if you need to.

By default only the alternate-allele calls are reported -- "which variants does
this subject carry". --hom-ref switches to every interrogated site -- "show me all
the sites for this subject" -- adding the reference calls to the same stream. The
gt column tells the two apart. --format vcf implies it, since a genotype matrix
that cannot tell 0/0 from ./. asserts far less than the data supports.

A 0/0 row means the whole genotype was reference. At a multiallelic record a 0/2
sample is not a carrier of allele 1, yet it is not reference either, so it
appears under neither -- writing 0/0 for it would be a genotype the source never
contained.

With no site restriction --hom-ref walks the entire sites catalog, since being
reference is a statement about every site the callset interrogated -- expect
roughly one row per variant in the source. --variant bounds it.

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

  --sample NAME       subject to report on (repeatable)
  --variant TARGET    a locus, region, contig, or a file of them (repeatable)
  --min-dp N          minimum depth for a call to count
  --hom-ref           report every interrogated site, not only the ALT calls
  --format F          tsv (default), json, vcf, or list. vcf emits a genotype
                      matrix, one record per site and one column per sample,
                      with AC/AN/AF/NS/nhomalt recomputed over those samples
  --store KIND        force the backend: vcf or parquet
  -v, --verbose       report the backend, the store's conversion settings, and
                      whether the quality gate could actually act, on stderr

Every non-VCF output uses one layout, whichever axes were named, so results can
be concatenated, sorted and cut the same way:

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
		// At least one selector. An empty Query is legal in the library -- it means
		// the whole store -- but making it the accidental result of a typo at the
		// command line is not worth the convenience.
		if len(vcfVarQuerySamples) == 0 && len(vcfVarQueryVariants) == 0 {
			return fmt.Errorf("give at least one --sample or --variant")
		}
		format := strings.ToLower(vcfVarQueryFormat)
		switch format {
		case "tsv", "json", "vcf":
		case "list":
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

		q, targets, err := buildQuery()
		if err != nil {
			return err
		}
		if vcfVarQueryVerbose {
			targets.report(cmd.ErrOrStderr())
		}
		// A genotype matrix has to distinguish an observed reference call from an
		// unobserved sample, or every non-carrier becomes ./. and the file asserts
		// far less than the data supports. So --format vcf implies --hom-ref, and
		// inherits its refusal on a store that cannot reconstruct reference calls.
		if format == "vcf" {
			q.IncludeRef = true
		}

		store, err := openVarStore(args[0], vcfVarQueryStore)
		if err != nil {
			return err
		}
		defer store.Close()

		if vcfVarQueryVerbose {
			describeStore(cmd, store, args[0], q.Gate)
		}
		warnUnknownSites(cmd, store, q)

		w, closeFn, err := openOutput(cmd, vcfVarQueryOutput)
		if err != nil {
			return err
		}
		out := bufio.NewWriter(w)

		var t *tally
		switch format {
		case "json":
			t, err = writeCallsJSON(out, store, q)
		case "vcf":
			t, err = writeCallsVCF(out, store, q, args[0])
		case "list":
			t, err = writeCallsList(out, store, q)
		default:
			t, err = writeCallsTSV(out, store, q, args[0])
		}
		if err != nil {
			return err
		}
		if vcfVarQueryVerbose {
			reportQuery(cmd, store, q, t)
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

// warnUnknownSites notes any queried variant the source never reported.
//
// Without this, such a variant returns zero carriers -- or all not-assayed
// under --hom-ref -- which reads exactly like a real negative result. The
// source simply never looked there, and only a gVCF could say otherwise.
func warnUnknownSites(cmd *cobra.Command, store varstore.Store, q varstore.Query) {
	for _, locus := range q.Loci {
		known, err := store.SiteKnown(locus)
		if err != nil || known {
			continue
		}
		fmt.Fprintf(cmd.ErrOrStderr(),
			"warning: %s is not in the source; reporting not-assayed for every sample.\n"+
				"         A VCF only supports queries for the variants it contains.\n", locus)
	}
}

// buildQuery turns the flags into one Query.
//
// Site selection and sample selection are independent axes, so any combination is
// legal: naming both is the variants-by-samples question that used to be a hard
// error, naming only --region asks for a window, and naming only --sample asks
// what that subject carries.
func buildQuery() (varstore.Query, *targetSet, error) {
	q := varstore.Query{
		Samples:    vcfVarQuerySamples,
		Gate:       varstore.Gate{MinDP: int32(vcfVarQueryMinDP)},
		IncludeRef: vcfVarQueryHomRef,
	}
	t, err := parseTargets(vcfVarQueryVariants)
	if err != nil {
		return q, nil, err
	}
	q.Loci, q.Spans = t.Loci, t.Spans
	return q, t, nil
}

// tally counts what a query emitted, for the verbose report. Gathered while
// streaming, since a streaming query cannot be counted afterwards.
type tally struct {
	alt, ref int
	withDP   int                       // rows carrying a real DP
	byLocus  map[varstore.Locus][2]int // [alt, ref] per locus
}

func newTally() *tally { return &tally{byLocus: map[varstore.Locus][2]int{}} }

func (t *tally) add(c varstore.Call) {
	if c.DP != varstore.Missing {
		t.withDP++
	}
	n := t.byLocus[c.Locus()]
	if varstore.IsAltCarrier(c.GT) {
		t.alt++
		n[0]++
	} else {
		t.ref++
		n[1]++
	}
	t.byLocus[c.Locus()] = n
}

// streamCalls runs the query and hands each row to fn, tallying as it goes.
func streamCalls(store varstore.Store, q varstore.Query, fn func(varstore.Call) error) (*tally, error) {
	seq, err := store.Calls(q)
	if err != nil {
		return nil, err
	}
	t := newTally()
	for c, err := range seq {
		if err != nil {
			return t, err
		}
		t.add(c)
		if err := fn(c); err != nil {
			return t, err
		}
	}
	return t, nil
}

// writeCallsTSV streams the shared tabular layout.
func writeCallsTSV(out *bufio.Writer, store varstore.Store, q varstore.Query, source string) (*tally, error) {
	provenance(out, source)
	fmt.Fprintln(out, strings.Join(tsvColumns, "\t"))
	minDP := vouchedMinDP(store)
	return streamCalls(store, q, func(c varstore.Call) error {
		_, err := fmt.Fprintf(out, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
			c.Chrom, c.Pos, c.Ref, c.Alt, c.SampleID, c.GT,
			intOrDot(c.DP), minDPOf(c, minDP),
			intOrDot(c.ADRef), intOrDot(c.ADAlt), intOrDot(c.GQ))
		return err
	})
}

// writeCallsJSON streams a JSON array rather than buffering one, so the shape is
// the same at any scale.
func writeCallsJSON(out *bufio.Writer, store varstore.Store, q varstore.Query) (*tally, error) {
	fmt.Fprintln(out, "[")
	first := true
	t, err := streamCalls(store, q, func(c varstore.Call) error {
		if !first {
			fmt.Fprintln(out, ",")
		}
		first = false
		b, err := json.Marshal(c)
		if err != nil {
			return err
		}
		_, err = out.Write(append([]byte("  "), b...))
		return err
	})
	if err != nil {
		return t, err
	}
	if !first {
		fmt.Fprintln(out)
	}
	fmt.Fprintln(out, "]")
	return t, nil
}

// writeCallsVCF streams a genotype matrix: one record per locus, one column per
// sample.
//
// It buffers only ONE locus at a time. Rows arrive in the store's order, so every
// row for a locus is contiguous -- the buffer is bounded by the sample count, not
// by the query, and the records come out in contig order without a comparator.
// The previous implementation collected every locus into a map and sorted at the
// end, which both bounded the output by memory and sorted chromosomes
// lexicographically (chr10 before chr2), producing a file tabix cannot index.
func writeCallsVCF(out *bufio.Writer, store varstore.Store, q varstore.Query, source string) (*tally, error) {
	names := q.Samples
	if len(names) == 0 {
		roster, err := store.Samples()
		if err != nil {
			return nil, err
		}
		names = roster
	}
	col := make(map[string]int, len(names))
	for i, n := range names {
		col[n] = i
	}

	fmt.Fprintln(out, "##fileformat=VCFv4.2")
	fmt.Fprintln(out, "##source="+buildinfo.String())
	fmt.Fprintln(out, "##cgkit_vcf-varqueryCommand="+buildinfo.CommandLine())
	fmt.Fprintln(out, "##cgkit_vcf-varquerySource="+source)
	fmt.Fprintln(out, `##INFO=<ID=AC,Number=A,Type=Integer,Description="Alt alleles among the samples in this file">`)
	fmt.Fprintln(out, `##INFO=<ID=AN,Number=1,Type=Integer,Description="Called alleles among the samples in this file (depth-gated)">`)
	fmt.Fprintln(out, `##INFO=<ID=AF,Number=A,Type=Float,Description="AC/AN">`)
	fmt.Fprintln(out, `##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with a call">`)
	fmt.Fprintln(out, `##INFO=<ID=nhomalt,Number=A,Type=Integer,Description="Samples homozygous for the alt allele">`)
	fmt.Fprintln(out, `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`)
	fmt.Fprintln(out, `##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">`)
	fmt.Fprintln(out, `##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">`)
	// AC/AN are recomputed over the samples in THIS file, not copied from the
	// store, so they stay correct over a sample subset. AN is depth-gated, because
	// what a store can vouch for is the samples inside its callable runs -- so it
	// can differ from the store's own ungated AN wherever a sample was called
	// below the threshold.
	fmt.Fprintln(out, "##cgkit_vcf-varqueryNote=AC/AN/AF/NS/nhomalt are recomputed over the samples in "+
		"this file; AN counts only samples the source can vouch for as called at the depth threshold")
	fmt.Fprintln(out, "#"+strings.Join(append([]string{
		"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
	}, names...), "\t"))

	var cur varstore.Locus
	row := make([]varstore.Call, len(names))
	held := false

	flush := func() error {
		if !held {
			return nil
		}
		var ac, an, ns, nhom int
		for i := range row {
			c := row[i]
			if c.GT == "" {
				continue // this sample had no row here: unobserved
			}
			ns++
			alleles, alt := countAlleles(c.GT)
			an += alleles
			ac += alt
			if alt >= 2 {
				nhom++
			}
		}
		af := "."
		if an > 0 {
			af = fmt.Sprintf("%.6g", float64(ac)/float64(an))
		}
		info := fmt.Sprintf("AC=%d;AN=%d;AF=%s;NS=%d;nhomalt=%d", ac, an, af, ns, nhom)

		cols := []string{cur.Chrom, fmt.Sprint(cur.Pos), ".", cur.Ref, cur.Alt, ".", ".", info, "GT:DP:GQ"}
		for i := range row {
			if row[i].GT == "" {
				cols = append(cols, "./.:.:.")
				continue
			}
			cols = append(cols, fmt.Sprintf("%s:%s:%s",
				row[i].GT, intOrDot(row[i].DP), intOrDot(row[i].GQ)))
		}
		_, err := fmt.Fprintln(out, strings.Join(cols, "\t"))
		for i := range row {
			row[i] = varstore.Call{}
		}
		held = false
		return err
	}

	t, err := streamCalls(store, q, func(c varstore.Call) error {
		if held && !varstore.SameLocus(c.Locus(), cur) {
			if err := flush(); err != nil {
				return err
			}
		}
		cur, held = c.Locus(), true
		if i, ok := col[c.SampleID]; ok {
			row[i] = c
		}
		return nil
	})
	if err != nil {
		return t, err
	}
	return t, flush()
}

// countAlleles returns how many alleles of a genotype were called, and how many
// of those are the alternate. Missing alleles ('.') count as neither, which is
// what makes AN a count of what was actually observed.
func countAlleles(gt string) (called, alt int) {
	for _, a := range strings.Split(strings.ReplaceAll(gt, "|", "/"), "/") {
		switch a {
		case "", ".":
		case "0":
			called++
		default:
			called++
			alt++
		}
	}
	return called, alt
}

// writeCallsList streams the distinct sample ids that had any row.
func writeCallsList(out *bufio.Writer, store varstore.Store, q varstore.Query) (*tally, error) {
	seen := map[string]bool{}
	return streamCalls(store, q, func(c varstore.Call) error {
		if seen[c.SampleID] {
			return nil
		}
		seen[c.SampleID] = true
		_, err := fmt.Fprintln(out, c.SampleID)
		return err
	})
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
	f.IntVar(&vcfVarQueryMinDP, "min-dp", 0, "Minimum DP for a call to count")
	f.BoolVar(&vcfVarQueryHomRef, "hom-ref", false, "Also report reference (0/0) calls, not only alternate carriers")
	f.StringVar(&vcfVarQueryFormat, "format", "tsv", "Output format: tsv, json, vcf, or list")
	f.StringVar(&vcfVarQueryStore, "store", "", "Force the backend: vcf or parquet (default: infer from the path)")
	f.BoolVarP(&vcfVarQueryVerbose, "verbose", "v", false, "Report the backend, store provenance and gate effect on stderr")
}

// reportQuery summarises what the query returned, from counts gathered while
// streaming -- a streaming query cannot be counted after the fact.
//
// Named loci get a line each, since that is the interactive case. Anything else
// gets totals: at panel scale a paragraph per locus would BE the output.
func reportQuery(cmd *cobra.Command, store varstore.Store, q varstore.Query, t *tally) {
	out := cmd.ErrOrStderr()
	if t == nil {
		return
	}

	const perLocusLimit = 32
	if len(q.Loci) > 0 && len(q.Loci) <= perLocusLimit {
		for _, l := range q.Loci {
			fmt.Fprintf(out, "variant  %s\n", l)
			// Allele frequency from the catalog, which knows the site even when
			// nobody carries it -- so this stays meaningful at AC 0. Note these are
			// the store's counts over ALL its samples, not the query's subset.
			if ps, ok := store.(*varstore.ParquetStore); ok {
				if site, found, err := ps.Site(l); err == nil && found {
					fmt.Fprintf(out, "  site        AC=%d AN=%d AF=%.6g  n_carriers=%d n_called=%d n_lowdp=%d\n",
						site.AC, site.AN, site.AF(), site.NCarriers, site.NCalled, site.NLowDP)
				}
			}
			n := t.byLocus[varstore.Locus{Chrom: l.Chrom, Pos: l.Pos, Ref: l.Ref, Alt: l.Alt}]
			if n == [2]int{} {
				// The store echoes its own contig spelling, which may differ from
				// the one asked for, so fall back to matching canonically.
				for k, v := range t.byLocus {
					if varstore.SameLocus(k, l) {
						n = v
						break
					}
				}
			}
			fmt.Fprintf(out, "  alt calls   %d\n", n[0])
			if q.IncludeRef {
				fmt.Fprintf(out, "  hom-ref     %d\n", n[1])
			}
		}
	}

	fmt.Fprintf(out, "total    %d alt call(s)", t.alt)
	if q.IncludeRef {
		fmt.Fprintf(out, ", %d reference call(s)", t.ref)
	}
	fmt.Fprintf(out, " over %d site(s)\n", len(t.byLocus))
	reportGate(cmd, store, q, t)
}

// reportGate says whether the quality gate could actually act.
//
// This is the diagnostic most worth surfacing. A gate over a field the data lacks
// admits everything rather than rejecting it, so a --min-dp that looks like a
// filter can be doing nothing at all -- and the numbers look equally plausible
// either way. Only the field census distinguishes them.
//
// The count of what the gate excluded needs the same query run ungated, so it is
// only offered when the query named loci: that is the interactive case, and at
// panel scale a second full pass for a verbose line is not worth it.
func reportGate(cmd *cobra.Command, store varstore.Store, q varstore.Query, t *tally) {
	if q.Gate.IsZero() {
		return
	}
	out := cmd.ErrOrStderr()
	rows := t.alt + t.ref

	if len(q.Loci) > 0 && len(q.Loci) <= 32 {
		ungated := q
		ungated.Gate = varstore.Gate{}
		if all, err := varstore.CollectCalls(store, ungated); err == nil {
			if n := len(all) - rows; n > 0 {
				fmt.Fprintf(out, "gate     excluded %d call(s)\n", n)
			} else {
				fmt.Fprintf(out, "gate     excluded 0 call(s)\n")
			}
		}
	}
	if q.Gate.MinDP > 0 && t.withDP == 0 && rows > 0 {
		fmt.Fprintf(out, "  WARNING: --min-dp %d had no effect; no call here carries DP\n", q.Gate.MinDP)
	}
}
