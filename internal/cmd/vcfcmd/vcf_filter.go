package vcfcmd

import (
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/compgenlab/hts/vcf"
	"github.com/compgenlab/hts/vcf/filter"
	"github.com/spf13/cobra"
)

var (
	vcfFilterOutput  string
	vcfFilterPassing bool
	vcfFilterFailing bool
	vcfFilterChain   []chainArg
)

var vcfFilterCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-filter <input.vcf>",
	Short:       "Filter a VCF file by stamping FILTER codes",
	Long: `Filter a VCF file. Each flag adds a filter; every variant runs through the
filters in command-line order, and a variant that fails a filter gets that
filter's code stamped into its FILTER column. A variant with no codes is PASS.

  --chrom-pass CSV   flag variants not on these chromosomes
  --chrom-fail CSV   flag variants on these chromosomes
  --snv / --indel    flag SNVs / indels
  --max-ins N        flag insertions longer than N
  --max-del N        flag deletions longer than N
  --qual FLOAT       flag variants with QUAL below FLOAT
  --het / --hom      flag heterozygous / homozygous variants (requires GT)

By default every variant is written (with FILTER updated). --passing writes only
variants that pass all filters; --failing writes only variants that fail one.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		chain, err := buildFilterChain()
		if err != nil {
			return err
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
		if err := chain.SetupHeaders(header); err != nil {
			return err
		}
		stampVcfProvenance(header, "vcf-filter")

		var writer *vcf.VcfWriter
		var closeFile func() error
		if vcfFilterOutput == "" || vcfFilterOutput == "-" {
			writer = vcf.NewVcfWriter(cmd.OutOrStdout())
		} else {
			w, err := vcf.OpenVcfWriter(vcfFilterOutput)
			if err != nil {
				return err
			}
			writer = w
			closeFile = w.Close
		}
		if err := writer.WriteHeader(header); err != nil {
			return err
		}

		for {
			rec, err := reader.NextRecord()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if err := chain.Apply(rec); err != nil {
				return err
			}
			filtered := rec.IsFiltered()
			if vcfFilterPassing && filtered {
				continue
			}
			if vcfFilterFailing && !filtered {
				continue
			}
			if err := writer.WriteRecord(rec); err != nil {
				return err
			}
		}
		if err := chain.Close(); err != nil {
			return err
		}
		if closeFile != nil {
			return closeFile()
		}
		return writer.Close()
	},
}

// buildFilterChain assembles the filter chain from the flags in command-line
// order.
func buildFilterChain() (*filter.Chain, error) {
	c := filter.NewChain()
	for _, fl := range vcfFilterChain {
		switch fl.kind {
		case "snv":
			c.Add(filter.NewSNV())
		case "indel":
			c.Add(filter.NewIndel())
		case "het":
			c.Add(filter.NewHeterozygous())
		case "hom":
			c.Add(filter.NewHomozygous())
		case "chrom-pass":
			c.Add(filter.NewChromPass(strings.Split(fl.value, ",")))
		case "chrom-fail":
			c.Add(filter.NewChromFail(strings.Split(fl.value, ",")))
		case "max-ins":
			n, err := strconv.Atoi(fl.value)
			if err != nil {
				return nil, fmt.Errorf("--max-ins: %w", err)
			}
			c.Add(filter.NewMaxIns(n))
		case "max-del":
			n, err := strconv.Atoi(fl.value)
			if err != nil {
				return nil, fmt.Errorf("--max-del: %w", err)
			}
			c.Add(filter.NewMaxDel(n))
		case "qual":
			q, err := strconv.ParseFloat(fl.value, 64)
			if err != nil {
				return nil, fmt.Errorf("--qual: %w", err)
			}
			c.Add(filter.NewQual(q))
		}
	}
	return c, nil
}

func init() {
	f := vcfFilterCmd.Flags()
	f.StringVarP(&vcfFilterOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
	f.BoolVar(&vcfFilterPassing, "passing", false, "Only output passing variants")
	f.BoolVar(&vcfFilterFailing, "failing", false, "Only output failing variants")
	registerChainVal(f, &vcfFilterChain, "chrom-pass", "Flag variants not on these chromosomes (CSV)")
	registerChainVal(f, &vcfFilterChain, "chrom-fail", "Flag variants on these chromosomes (CSV)")
	registerChainBool(f, &vcfFilterChain, "snv", "Flag SNVs")
	registerChainBool(f, &vcfFilterChain, "indel", "Flag indels")
	registerChainVal(f, &vcfFilterChain, "max-ins", "Flag insertions longer than this length")
	registerChainVal(f, &vcfFilterChain, "max-del", "Flag deletions longer than this length")
	registerChainVal(f, &vcfFilterChain, "qual", "Flag variants with QUAL below this value")
	registerChainBool(f, &vcfFilterChain, "het", "Flag heterozygous variants (requires GT)")
	registerChainBool(f, &vcfFilterChain, "hom", "Flag homozygous variants (requires GT)")
}
