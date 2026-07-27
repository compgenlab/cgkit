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
  --store KIND        force the backend: vcf or parquet`,
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

		w, closeFn, err := openOutput(cmd, vcfVarQueryOutput)
		if err != nil {
			return err
		}
		out := bufio.NewWriter(w)

		if nSample > 0 {
			err = runVarQuerySamples(out, store, span, gate, format, args[0])
		} else {
			err = runVarQueryVariants(out, store, gate, format, args[0])
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
func runVarQueryVariants(out *bufio.Writer, store varstore.Store,
	gate varstore.Gate, format, source string) error {

	type variantResult struct {
		Variant string                 `json:"variant"`
		States  []varstore.SampleState `json:"states,omitempty"`
		Calls   []varstore.Call        `json:"calls,omitempty"`
	}
	var results []variantResult

	for _, v := range vcfVarQueryVariants {
		locus, err := varstore.ParseLocus(v)
		if err != nil {
			return err
		}
		r := variantResult{Variant: locus.String()}
		if vcfVarQueryClassify {
			states, err := store.Classify(locus, gate)
			if err != nil {
				return err
			}
			r.States = states
		} else {
			calls, err := store.Carriers(locus, gate)
			if err != nil {
				return err
			}
			varstore.SortCalls(calls)
			r.Calls = calls
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
	if vcfVarQueryClassify {
		fmt.Fprintln(out, strings.Join([]string{"variant", "sample", "state", "gt", "dp", "gq"}, "\t"))
		for _, r := range results {
			for _, s := range r.States {
				gt, dp, gq := ".", ".", "."
				if s.Call != nil {
					gt, dp, gq = s.Call.GT, intOrDot(s.Call.DP), intOrDot(s.Call.GQ)
				}
				fmt.Fprintf(out, "%s\t%s\t%s\t%s\t%s\t%s\n", r.Variant, s.SampleID, s.State, gt, dp, gq)
			}
		}
		return nil
	}
	fmt.Fprintln(out, strings.Join([]string{
		"variant", "sample", "gt", "dp", "ad_ref", "ad_alt", "gq",
	}, "\t"))
	for _, r := range results {
		for _, c := range r.Calls {
			fmt.Fprintf(out, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", r.Variant, c.SampleID, c.GT,
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
}
