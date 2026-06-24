package vcfcmd

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/compgenlab/hts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfExportOutput      string
	vcfExportID          bool
	vcfExportQual        bool
	vcfExportFilter      bool
	vcfExportInfo        []string
	vcfExportFormat      []string
	vcfExportNoHeader    bool
	vcfExportNoVCFHeader bool
	vcfExportPassing     bool
	vcfExportOnlySNVs    bool
	vcfExportOnlyIndels  bool
	vcfExportMissBlank   bool
	vcfExportRegion      string
)

// exporter produces one or more output columns per record.
type exporter interface {
	fieldNames() []string
	export(rec *vcf.VcfRecord, out *[]string) error
}

var vcfExportCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-export <input.vcf>",
	Short:       "Export information from a VCF file as a tab-delimited file",
	Long: `Export information from a VCF file as a tab-delimited file.

Each output row begins with chrom, pos, ref, alt, followed by the requested
columns. --id, --qual, and --filter export the VCF ID/QUAL/FILTER columns.
--info and --format export INFO and FORMAT fields; both accept glob patterns
(* and ?) for the field name and may be repeated or comma-separated.

For --info the form is KEY[:ALLELE]; for --format it is
ID[:SAMPLE[:ALLELE[:NEWSAMPLE]]]. ALLELE may be one of sum, min, max, ref, or
alt1 (or a numeric index), or left blank for the whole value.

Columns are emitted in this order: id, qual, filter, then each --info (in the
order given), then each --format.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if vcfExportOnlySNVs && vcfExportOnlyIndels {
			return fmt.Errorf("you can't set both --only-snvs and --only-indels")
		}

		src, err := openRecordSource(cmd, args[0], vcfExportRegion)
		if err != nil {
			return err
		}
		defer src.close()

		header := src.header

		missing := "."
		if vcfExportMissBlank {
			missing = ""
		}

		chain, err := buildExportChain(header, missing)
		if err != nil {
			return err
		}

		w, closeFn, err := openOutput(cmd, vcfExportOutput)
		if err != nil {
			return err
		}
		out := bufio.NewWriter(w)

		if !vcfExportNoVCFHeader {
			for _, line := range header.MetaLines() {
				fmt.Fprintln(out, line)
			}
		}

		if !vcfExportNoHeader {
			cols := []string{"chrom", "pos", "ref", "alt"}
			for _, e := range chain {
				cols = append(cols, e.fieldNames()...)
			}
			fmt.Fprintln(out, strings.Join(cols, "\t"))
		}

		for {
			rec, err := src.next()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if vcfExportPassing && rec.IsFiltered() {
				continue
			}
			if vcfExportOnlySNVs && rec.IsIndel() {
				continue
			}
			if vcfExportOnlyIndels && !rec.IsIndel() {
				continue
			}

			outs := []string{rec.Chrom, strconv.Itoa(rec.Pos), rec.Ref, strings.Join(rec.Alt(), ",")}
			for _, e := range chain {
				if err := e.export(rec, &outs); err != nil {
					return err
				}
			}
			fmt.Fprintln(out, strings.Join(outs, "\t"))
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

// buildExportChain assembles the exporters in a fixed column order: id, qual,
// filter, then each --info, then each --format.
func buildExportChain(header *vcf.VcfHeader, missing string) ([]exporter, error) {
	var chain []exporter
	if vcfExportID {
		chain = append(chain, &idExporter{missing: missing})
	}
	if vcfExportQual {
		chain = append(chain, &qualExporter{missing: missing})
	}
	if vcfExportFilter {
		chain = append(chain, &filterExporter{})
	}

	exported := map[[2]string]bool{}
	for _, raw := range vcfExportInfo {
		for _, val := range strings.Split(raw, ",") {
			key, allele, ignoreMissing := parseFieldSpec(val)
			ids := []string{}
			for _, id := range header.MatchInfoIDs(key) {
				k := [2]string{id, allele}
				if exported[k] {
					continue
				}
				exported[k] = true
				ids = append(ids, id)
			}
			chain = append(chain, &infoExporter{allele: allele, ignoreMissing: ignoreMissing, missing: missing, ids: ids})
		}
	}

	for _, raw := range vcfExportFormat {
		for _, val := range strings.Split(raw, ",") {
			spec, err := parseFormatSpec(val)
			if err != nil {
				return nil, err
			}
			spec.missing = missing
			spec.ids = header.MatchFormatIDs(spec.key)
			spec.samples = header.Samples()
			if spec.sample != "" {
				spec.sampleIdx = header.SampleIndex(spec.sample)
			} else {
				spec.sampleIdx = -1
			}
			chain = append(chain, spec)
		}
	}
	return chain, nil
}

// parseFieldSpec splits an --info value into key, allele, and ignoreMissing,
// honoring the trailing ":!" (require) / ":?" (ignore) markers.
func parseFieldSpec(val string) (key, allele string, ignoreMissing bool) {
	ignoreMissing = true
	if strings.HasSuffix(val, ":!") {
		ignoreMissing = false
		val = val[:len(val)-2]
	} else if strings.HasSuffix(val, ":?") {
		ignoreMissing = true
		val = val[:len(val)-2]
	}
	spl := strings.Split(val, ":")
	key = spl[0]
	if len(spl) > 1 {
		allele = spl[1]
	}
	return key, allele, ignoreMissing
}

func parseFormatSpec(val string) (*formatExporter, error) {
	ignoreMissing := true
	if strings.HasSuffix(val, ":!") {
		ignoreMissing = false
		val = val[:len(val)-2]
	} else if strings.HasSuffix(val, ":?") {
		ignoreMissing = true
		val = val[:len(val)-2]
	}
	spl := strings.Split(val, ":")
	var key, sample, allele, newName string
	if len(spl) > 0 {
		key = spl[0]
	}
	if len(spl) > 1 {
		sample = spl[1]
	}
	if len(spl) > 2 {
		allele = spl[2]
	}
	if len(spl) > 3 {
		newName = spl[3]
	}
	if key == "" {
		return nil, fmt.Errorf("missing field name for --format")
	}
	if newName != "" && sample == "" {
		return nil, fmt.Errorf("--format with a new sample name requires a SAMPLE")
	}
	return &formatExporter{key: key, sample: sample, allele: allele, ignoreMissing: ignoreMissing, newName: newName}, nil
}

type idExporter struct{ missing string }

func (e *idExporter) fieldNames() []string { return []string{"ID"} }
func (e *idExporter) export(rec *vcf.VcfRecord, out *[]string) error {
	id := rec.ID()
	if id == "" {
		*out = append(*out, e.missing)
	} else {
		*out = append(*out, id)
	}
	return nil
}

type qualExporter struct{ missing string }

func (e *qualExporter) fieldNames() []string { return []string{"QUAL"} }
func (e *qualExporter) export(rec *vcf.VcfRecord, out *[]string) error {
	q := rec.Qual()
	if q == -1 {
		*out = append(*out, e.missing)
	} else {
		*out = append(*out, javaDouble(q))
	}
	return nil
}

type filterExporter struct{}

func (e *filterExporter) fieldNames() []string { return []string{"FILTER"} }
func (e *filterExporter) export(rec *vcf.VcfRecord, out *[]string) error {
	if !rec.IsFiltered() {
		*out = append(*out, "PASS")
	} else {
		*out = append(*out, strings.Join(rec.Filters(), ","))
	}
	return nil
}

type infoExporter struct {
	allele        string
	ignoreMissing bool
	missing       string
	ids           []string
}

func (e *infoExporter) fieldNames() []string { return e.ids }
func (e *infoExporter) export(rec *vcf.VcfRecord, out *[]string) error {
	for _, id := range e.ids {
		v, ok := rec.Info().Get(id)
		if !ok {
			if e.ignoreMissing {
				*out = append(*out, "")
				continue
			}
			return fmt.Errorf("unable to find INFO field: %s", id)
		}
		if e.ignoreMissing && v.IsEmpty() {
			// A present flag exports as its key name.
			*out = append(*out, id)
			continue
		}
		s, emit, err := selectValue(v, e.allele, e.missing)
		if err != nil {
			return err
		}
		if emit {
			*out = append(*out, s)
		}
	}
	return nil
}

type formatExporter struct {
	key           string
	sample        string
	allele        string
	ignoreMissing bool
	missing       string
	newName       string

	ids       []string
	samples   []string
	sampleIdx int
}

func (e *formatExporter) fieldNames() []string {
	var out []string
	if e.sample != "" {
		label := e.sample
		if e.newName != "" {
			label = e.newName
		}
		for _, id := range e.ids {
			out = append(out, label+":"+id)
		}
	} else {
		for _, s := range e.samples {
			for _, id := range e.ids {
				out = append(out, s+":"+id)
			}
		}
	}
	return out
}

func (e *formatExporter) export(rec *vcf.VcfRecord, out *[]string) error {
	if rec.NumSamples() == 0 {
		return nil
	}
	if e.sampleIdx == -1 {
		for i := 0; i < rec.NumSamples(); i++ {
			attr, err := rec.Sample(i)
			if err != nil {
				return err
			}
			if err := e.exportVal(attr, out); err != nil {
				return err
			}
		}
		return nil
	}
	attr, err := rec.Sample(e.sampleIdx)
	if err != nil {
		return err
	}
	return e.exportVal(attr, out)
}

func (e *formatExporter) exportVal(attr *vcf.Attributes, out *[]string) error {
	for _, id := range e.ids {
		v, ok := attr.Get(id)
		if !ok {
			if e.ignoreMissing {
				*out = append(*out, "")
				continue
			}
			return fmt.Errorf("unable to find FORMAT field: %s", id)
		}
		s, emit, err := selectValue(v, e.allele, e.missing)
		if err != nil {
			return err
		}
		if emit {
			*out = append(*out, s)
		}
	}
	return nil
}

// selectValue renders an attribute value for the given allele selector. The
// boolean is false when nothing should be emitted (matching ngsutilsj, which
// drops a column when min/max cannot be computed).
func selectValue(v vcf.AttrValue, allele, missing string) (string, bool, error) {
	switch allele {
	case "sum":
		d, err := v.FloatFor("sum")
		if err != nil {
			return "", false, err
		}
		return formatStripped(d), true, nil
	case "min", "max":
		d, err := v.FloatFor(allele)
		if err != nil {
			return "", false, err
		}
		if math.IsNaN(d) {
			return "", false, nil
		}
		return formatStripped(d), true, nil
	default:
		if v.IsMissing() {
			return missing, true, nil
		}
		s, err := v.StringFor(allele)
		if err != nil {
			return "", false, err
		}
		return s, true, nil
	}
}

// javaDouble formats a float the way Java's Double.toString does for ordinary
// magnitudes: an integral value keeps a ".0" suffix (e.g. 50 -> "50.0").
func javaDouble(f float64) string {
	s := strconv.FormatFloat(f, 'f', -1, 64)
	if !strings.Contains(s, ".") {
		s += ".0"
	}
	return s
}

// formatStripped formats a float as a plain decimal, dropping a trailing ".0".
func formatStripped(f float64) string {
	s := strconv.FormatFloat(f, 'f', -1, 64)
	return strings.TrimSuffix(s, ".0")
}

func init() {
	vcfExportCmd.Flags().StringVarP(&vcfExportOutput, "output", "o", "-", "Output filename (- for stdout)")
	vcfExportCmd.Flags().BoolVar(&vcfExportID, "id", false, "Export the VCF ID column")
	vcfExportCmd.Flags().BoolVar(&vcfExportQual, "qual", false, "Export the VCF QUAL column")
	vcfExportCmd.Flags().BoolVar(&vcfExportFilter, "filter", false, "Export the VCF FILTER column")
	vcfExportCmd.Flags().StringArrayVar(&vcfExportInfo, "info", nil, "Export an INFO field: KEY[:ALLELE] (repeatable, glob, comma-separated)")
	vcfExportCmd.Flags().StringArrayVar(&vcfExportFormat, "format", nil, "Export a FORMAT field: ID[:SAMPLE[:ALLELE[:NEWSAMPLE]]] (repeatable, glob, comma-separated)")
	vcfExportCmd.Flags().BoolVar(&vcfExportNoHeader, "no-header", false, "Don't write the column-name row")
	vcfExportCmd.Flags().BoolVar(&vcfExportNoVCFHeader, "no-vcf-header", false, "Don't write the VCF metadata header lines")
	vcfExportCmd.Flags().BoolVar(&vcfExportPassing, "passing", false, "Only export passing variants")
	vcfExportCmd.Flags().BoolVar(&vcfExportOnlySNVs, "only-snvs", false, "Only export SNVs")
	vcfExportCmd.Flags().BoolVar(&vcfExportOnlyIndels, "only-indels", false, "Only export indels")
	vcfExportCmd.Flags().BoolVar(&vcfExportMissBlank, "missing-blank", false, "Render missing values as an empty string instead of \".\"")
	vcfExportCmd.Flags().StringVar(&vcfExportRegion, "region", "", "Only variants in this 1-based region (chrom:start-end, or chrom); requires a tabix-indexed file")
}
