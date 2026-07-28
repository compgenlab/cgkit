package vcfcmd

import (
	"bufio"
	"encoding/json"
	"fmt"
	"os"
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
	vcfVarQueryMinGQ    int
	vcfVarQueryClassify bool
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
a VCF (plain or bgzipped) or a Parquet store written by vcf-toparquet, given
either as its base name or as any one of its three files. The backend is
inferred from the path; override with --store.

Two modes, one of which must be given:

  --sample NAME     report the variants that subject carries
  --variant LOCUS   report the subjects carrying that variant (chrom:pos:ref:alt)

Both are repeatable. --min-dp and --min-gq gate calls by quality; a gate over
data lacking that field admits everything rather than rejecting it, since an
absent quality score is not evidence of a bad one.

With --classify, every sample is resolved to one of four states rather than
only the carriers being listed:

  carrier       an alternate call passing the gate
  uncertain     an alternate call below the gate
  non-carrier   observed at this position, and observed to be reference
  not-assayed   no observation to draw on

A VCF answers this directly, because it stores an explicit genotype for every
sample at every record. A Parquet store reconstructs it from its sites and
callable-regions files, and will refuse rather than guess if either is missing
or was built without coverage information -- reporting an unobserved sample as
a non-carrier would invent an observation.

  --sample NAME       subject to report variants for (repeatable)
  --variant LOCUS     variant to report carriers of (repeatable)
  --region R          restrict to a 1-based region (chrom:start-end)
  --min-dp N          minimum depth for a call to count
  --min-gq N          minimum genotype quality for a call to count
  --classify          resolve every sample, not just carriers
  --format F          tsv (default), json, vcf, or list
  --store KIND        force the backend: vcf or parquet
  -v, --verbose       report the backend, the store's conversion settings, and
                      whether the quality gate could actually act, on stderr

Tabular output splits the locus into four leading columns (chrom, pos, ref,
alt) rather than one packed chrom:pos:ref:alt field, so it can be cut and
sorted on position directly. The chromosome is echoed the way it was asked for,
whichever naming convention the file itself uses.`,
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
		default:
			return fmt.Errorf("unknown format %q (use tsv, json, vcf, or list)", format)
		}
		if vcfVarQueryClassify && nVariant == 0 {
			return fmt.Errorf("--classify applies to --variant mode, which resolves samples at a locus")
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
		gate := varstore.Gate{MinDP: int32(vcfVarQueryMinDP), MinGQ: int32(vcfVarQueryMinGQ)}
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
			fmt.Fprintf(out, "  callable    not tracked (--no-callable): --classify will refuse\n")
		} else {
			fmt.Fprintf(out, "  min-dp      %d at conversion\n", p.MinDP)
			if g.MinDP > 0 && g.MinDP != p.MinDP {
				fmt.Fprintf(out, "  NOTE: querying at --min-dp %d but the runs were built at %d;\n"+
					"        non-carrier vs not-assayed will not match a direct VCF query.\n",
					g.MinDP, p.MinDP)
			}
		}
	case *varstore.VcfStore:
		names, _ := s.Samples()
		fmt.Fprintf(out, "store    vcf %s (%d samples)\n", path, len(names))
	}
}

// reportGate says whether the quality gate could actually act.
//
// A gate over a field the data lacks admits everything rather than rejecting
// it, so a --min-gq that looks like a filter can be doing nothing at all. This
// is the one diagnostic most worth surfacing: the numbers look plausible either
// way, and only the field census distinguishes them.
func reportGate(cmd *cobra.Command, g varstore.Gate, calls []varstore.Call, excluded int) {
	if g.IsZero() {
		return
	}
	out := cmd.ErrOrStderr()
	var haveDP, haveGQ int
	for _, c := range calls {
		if c.DP != varstore.Missing {
			haveDP++
		}
		if c.GQ != varstore.Missing {
			haveGQ++
		}
	}
	fmt.Fprintf(out, "gate     excluded %d call(s)\n", excluded)
	if g.MinDP > 0 && haveDP == 0 && len(calls) > 0 {
		fmt.Fprintf(out, "  WARNING: --min-dp %d had no effect; no call here carries DP\n", g.MinDP)
	}
	if g.MinGQ > 0 && haveGQ == 0 && len(calls) > 0 {
		fmt.Fprintf(out, "  WARNING: --min-gq %d had no effect; no call here carries GQ\n", g.MinGQ)
	}
}

// warnUnknownSites notes any queried variant the source never reported.
//
// Without this, such a variant returns zero carriers -- or all not-assayed
// under --classify -- which reads exactly like a real negative result. The
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
		calls, err := store.Variants(s, span, gate)
		if err != nil {
			return err
		}
		varstore.SortCalls(calls)
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
	fmt.Fprintln(out, strings.Join([]string{
		"sample", "chrom", "pos", "ref", "alt", "gt", "dp", "ad_ref", "ad_alt", "gq",
	}, "\t"))
	for _, r := range results {
		for _, c := range r.Calls {
			fmt.Fprintf(out, "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				c.SampleID, c.Chrom, c.Pos, c.Ref, c.Alt, c.GT,
				intOrDot(c.DP), intOrDot(c.ADRef), intOrDot(c.ADAlt), intOrDot(c.GQ))
		}
	}
	return nil
}

// writeVarQueryVcf emits a minimal VCF of the carried variants. It is
// deliberately spare: the store keeps only alternate calls, so anything not
// carried by the requested samples cannot be reconstructed and is absent
// rather than written as reference.
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
		Variant string                 `json:"variant"`
		Locus   varstore.Locus         `json:"locus"`
		States  []varstore.SampleState `json:"states,omitempty"`
		Calls   []varstore.Call        `json:"calls,omitempty"`
	}
	var results []variantResult

	for _, v := range vcfVarQueryVariants {
		locus, err := varstore.ParseLocus(v)
		if err != nil {
			return err
		}
		r := variantResult{Variant: locus.String(), Locus: locus}

		// Compare against the ungated result so the gate's effect is a measured
		// number rather than an inference from the row count.
		var ungated []varstore.Call
		if vcfVarQueryVerbose && !gate.IsZero() {
			ungated, err = store.Carriers(locus, varstore.Gate{})
			if err != nil {
				return err
			}
		}

		if vcfVarQueryClassify {
			states, err := store.Classify(locus, gate)
			if err != nil {
				return err
			}
			r.States = states
			if vcfVarQueryVerbose {
				reportVariant(cmd, store, locus, len(states), tallyStates(states))
			}
		} else {
			calls, err := store.Carriers(locus, gate)
			if err != nil {
				return err
			}
			varstore.SortCalls(calls)
			r.Calls = calls
			if vcfVarQueryVerbose {
				reportVariant(cmd, store, locus, len(calls), nil)
			}
		}
		if vcfVarQueryVerbose && !gate.IsZero() {
			kept := len(r.Calls)
			if r.States != nil {
				kept = tallyStates(r.States)[varstore.StateCarrier]
			}
			reportGate(cmd, gate, ungated, len(ungated)-kept)
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
			for _, s := range r.States {
				if s.State == varstore.StateCarrier && !seen[s.SampleID] {
					seen[s.SampleID] = true
					fmt.Fprintln(out, s.SampleID)
				}
			}
		}
		return nil
	}

	provenance(out, source)
	// The locus is emitted as four columns rather than one packed
	// chrom:pos:ref:alt field, matching --sample mode and keeping the output
	// sortable and cut-able on position without re-splitting a composite key.
	if vcfVarQueryClassify {
		fmt.Fprintln(out, strings.Join([]string{
			"chrom", "pos", "ref", "alt", "sample", "state", "gt", "dp", "gq",
		}, "\t"))
		for _, r := range results {
			for _, s := range r.States {
				gt, dp, gq := ".", ".", "."
				if s.Call != nil {
					gt, dp, gq = s.Call.GT, intOrDot(s.Call.DP), intOrDot(s.Call.GQ)
				}
				fmt.Fprintf(out, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					r.Locus.Chrom, r.Locus.Pos, r.Locus.Ref, r.Locus.Alt,
					s.SampleID, s.State, gt, dp, gq)
			}
		}
		return nil
	}
	fmt.Fprintln(out, strings.Join([]string{
		"chrom", "pos", "ref", "alt", "sample", "gt", "dp", "ad_ref", "ad_alt", "gq",
	}, "\t"))
	for _, r := range results {
		for _, c := range r.Calls {
			// The queried locus, not the call's stored spelling: every row of a
			// --variant query then reads back the way it was asked for, and both
			// sub-modes agree even though only this one has calls to draw on.
			fmt.Fprintf(out, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
				r.Locus.Chrom, r.Locus.Pos, r.Locus.Ref, r.Locus.Alt,
				c.SampleID, c.GT,
				intOrDot(c.DP), intOrDot(c.ADRef), intOrDot(c.ADAlt), intOrDot(c.GQ))
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

func init() {
	f := vcfVarQueryCmd.Flags()
	f.StringVarP(&vcfVarQueryOutput, "output", "o", "-", "Output filename (- for stdout)")
	f.StringArrayVar(&vcfVarQuerySamples, "sample", nil, "Report variants carried by this subject (repeatable)")
	f.StringArrayVar(&vcfVarQueryVariants, "variant", nil, "Report subjects carrying this variant, as chrom:pos:ref:alt (repeatable)")
	f.StringVar(&vcfVarQueryRegion, "region", "", "Only variants in this 1-based region (chrom:start-end, or chrom); requires a tabix-indexed file")
	f.IntVar(&vcfVarQueryMinDP, "min-dp", 0, "Minimum DP for a call to count")
	f.IntVar(&vcfVarQueryMinGQ, "min-gq", 0, "Minimum GQ for a call to count")
	f.BoolVar(&vcfVarQueryClassify, "classify", false, "Resolve every sample to carrier/uncertain/non-carrier/not-assayed")
	f.StringVar(&vcfVarQueryFormat, "format", "tsv", "Output format: tsv, json, vcf, or list")
	f.StringVar(&vcfVarQueryStore, "store", "", "Force the backend: vcf or parquet (default: infer from the path)")
	f.BoolVarP(&vcfVarQueryVerbose, "verbose", "v", false, "Report the backend, store provenance and gate effect on stderr")
}

// tallyStates counts how many samples landed in each classification state.
func tallyStates(states []varstore.SampleState) map[varstore.State]int {
	out := map[varstore.State]int{}
	for _, s := range states {
		out[s.State]++
	}
	return out
}

// reportVariant summarises one queried locus: what the catalog knows about it
// and, under --classify, how the samples fell out.
func reportVariant(cmd *cobra.Command, store varstore.Store, l varstore.Locus,
	n int, tally map[varstore.State]int) {

	out := cmd.ErrOrStderr()
	fmt.Fprintf(cmd.ErrOrStderr(), "variant  %s\n", l)

	// Allele frequency comes from the catalog, which knows the site even when
	// nobody carries it -- so this stays meaningful at AC 0.
	if ps, ok := store.(*varstore.ParquetStore); ok {
		if site, found, err := ps.Site(l); err == nil && found {
			fmt.Fprintf(out, "  site        AC=%d AN=%d AF=%.6g  n_carriers=%d n_called=%d n_lowdp=%d\n",
				site.AC, site.AN, site.AF(), site.NCarriers, site.NCalled, site.NLowDP)
		}
	}
	if tally == nil {
		fmt.Fprintf(out, "  carriers    %d\n", n)
		return
	}
	fmt.Fprintf(out, "  classified  %d samples:", n)
	for _, st := range []varstore.State{
		varstore.StateCarrier, varstore.StateUncertain,
		varstore.StateNonCarrier, varstore.StateNotAssayed,
	} {
		fmt.Fprintf(out, "  %s=%d", st, tally[st])
	}
	fmt.Fprintln(out)
}
