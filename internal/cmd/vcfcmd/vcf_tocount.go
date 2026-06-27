package vcfcmd

import (
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/compgenlab/hts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfToCountOutput  string
	vcfToCountSample  string
	vcfToCountROAO    bool
	vcfToCountAF      bool
	vcfToCountTotal   bool
	vcfToCountHet     bool
	vcfToCountPassing bool
)

var vcfToCountCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-tocount <input.vcf>",
	Short:       "Convert a VCF to a count file using the AD (or RO/AO) format field",
	Long: `Write per-allele reference and alternate counts as a tab-delimited file, one
row per ALT allele. Counts come from the AD FORMAT field (or RO/AO with
--use-ro-ao).

  --sample NAME   sample to use (default: the first sample)
  --use-ro-ao     use RO/AO instead of AD
  --af            add an alt_freq column
  --total         add a total_count column
  --het           only count heterozygous variants (GT 0/1)
  --passing       only count passing variants`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		reader, err := openVcfInput(cmd, args[0])
		if err != nil {
			return err
		}
		defer reader.Close()

		header, err := reader.Header()
		if err != nil {
			return err
		}
		if err := requireCountFields(header); err != nil {
			return err
		}
		sampleIdx := 0
		if vcfToCountSample != "" {
			sampleIdx = header.SampleIndex(vcfToCountSample)
			if sampleIdx < 0 {
				return fmt.Errorf("sample not found: %s", vcfToCountSample)
			}
		}

		out, closeFn, err := openOutput(cmd, vcfToCountOutput)
		if err != nil {
			return err
		}
		fmt.Fprintln(out, "## program: "+buildinfo.String())
		fmt.Fprintln(out, "## cmd: "+buildinfo.CommandLine())
		fmt.Fprintln(out, "## vcf-input: "+args[0])
		hdr := []string{"chrom", "pos", "ref", "alt", "ref_count", "alt_count"}
		if vcfToCountAF {
			hdr = append(hdr, "alt_freq")
		}
		if vcfToCountTotal {
			hdr = append(hdr, "total_count")
		}
		fmt.Fprintln(out, strings.Join(hdr, "\t"))

		for {
			rec, err := reader.NextRecord()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if vcfToCountPassing && rec.IsFiltered() {
				continue
			}
			if err := writeCounts(out, rec, sampleIdx); err != nil {
				return err
			}
		}
		if closeFn != nil {
			return closeFn()
		}
		return nil
	},
}

func requireCountFields(header *vcf.VcfHeader) error {
	has := func(id string) bool { _, ok := header.FormatDef(id); return ok }
	if !vcfToCountROAO {
		if !has("AD") {
			return fmt.Errorf("the VCF must contain the AD format field")
		}
		return nil
	}
	if !has("RO") || !has("AO") {
		return fmt.Errorf("the VCF must contain the RO and AO format fields")
	}
	return nil
}

func writeCounts(out io.Writer, rec *vcf.VcfRecord, sampleIdx int) error {
	sample, err := rec.Sample(sampleIdx)
	if err != nil {
		return err
	}
	if vcfToCountHet {
		gt, ok := sample.Get("GT")
		if !ok {
			return fmt.Errorf("missing GT field (%s:%d)", rec.Chrom, rec.Pos)
		}
		if v := gt.String(); v != "0/1" && v != "0|1" {
			return nil
		}
	}
	for i, alt := range rec.Alt() {
		refCount, altCount, total, err := alleleCounts(rec, sample, i)
		if err != nil {
			return err
		}
		row := []string{rec.Chrom, strconv.Itoa(rec.Pos), rec.Ref, alt,
			strconv.Itoa(refCount), strconv.Itoa(altCount)}
		if vcfToCountAF {
			row = append(row, javaDouble(float64(altCount)/float64(altCount+refCount)))
		}
		if vcfToCountTotal {
			row = append(row, strconv.Itoa(total))
		}
		fmt.Fprintln(out, strings.Join(row, "\t"))
	}
	return nil
}

func alleleCounts(rec *vcf.VcfRecord, sample *vcf.Attributes, i int) (refCount, altCount, total int, err error) {
	if !vcfToCountROAO {
		ad, ok := sample.Get("AD")
		if !ok {
			return 0, 0, 0, fmt.Errorf("missing AD field (%s:%d)", rec.Chrom, rec.Pos)
		}
		vals, err := parseInts(strings.Split(ad.String(), ","))
		if err != nil || len(vals) < i+2 {
			return 0, 0, 0, fmt.Errorf("invalid AD field %q (%s:%d)", ad.String(), rec.Chrom, rec.Pos)
		}
		refCount, altCount = vals[0], vals[i+1]
		for _, v := range vals {
			total += v
		}
		return refCount, altCount, total, nil
	}
	ro, ok := sample.Get("RO")
	if !ok {
		return 0, 0, 0, fmt.Errorf("missing RO field (%s:%d)", rec.Chrom, rec.Pos)
	}
	ao, ok := sample.Get("AO")
	if !ok {
		return 0, 0, 0, fmt.Errorf("missing AO field (%s:%d)", rec.Chrom, rec.Pos)
	}
	refCount, err = strconv.Atoi(ro.String())
	if err != nil {
		return 0, 0, 0, fmt.Errorf("invalid RO field %q (%s:%d)", ro.String(), rec.Chrom, rec.Pos)
	}
	aoVals, err := parseInts(strings.Split(ao.String(), ","))
	if err != nil || len(aoVals) <= i {
		return 0, 0, 0, fmt.Errorf("invalid AO field %q (%s:%d)", ao.String(), rec.Chrom, rec.Pos)
	}
	altCount = aoVals[i]
	total = refCount
	for _, v := range aoVals {
		total += v
	}
	return refCount, altCount, total, nil
}

func parseInts(fields []string) ([]int, error) {
	out := make([]int, len(fields))
	for i, f := range fields {
		n, err := strconv.Atoi(f)
		if err != nil {
			return nil, err
		}
		out[i] = n
	}
	return out, nil
}

func init() {
	f := vcfToCountCmd.Flags()
	f.StringVarP(&vcfToCountOutput, "output", "o", "-", "Output filename (- for stdout)")
	f.StringVar(&vcfToCountSample, "sample", "", "Sample to use (default: the first sample)")
	f.BoolVar(&vcfToCountROAO, "use-ro-ao", false, "Use RO/AO format fields instead of AD")
	f.BoolVar(&vcfToCountAF, "af", false, "Output alternate allele frequency")
	f.BoolVar(&vcfToCountTotal, "total", false, "Output total allele count")
	f.BoolVar(&vcfToCountHet, "het", false, "Only count heterozygous variants (GT 0/1)")
	f.BoolVar(&vcfToCountPassing, "passing", false, "Only count passing variants")
}
