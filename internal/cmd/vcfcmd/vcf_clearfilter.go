package vcfcmd

import (
	"io"
	"strings"

	"github.com/compgenlab/hts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfClearFilterOutput  string
	vcfClearFilterFilters []string
	vcfClearFilterOnly    bool
	vcfClearFilterPassing bool
)

var vcfClearFilterCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-clearfilter <input.vcf>",
	Short:       "Remove a filter from a VCF file",
	Long: `Remove previously set FILTER codes from variants.

By default this removes the named filter(s) from every variant. With --only, a
filter is only cleared when the named filters are the *only* codes on a variant.
Cleared codes are recorded in the CG_CLEARED_FILTER INFO field.

  --filter VAL   filter code to clear (comma-separated, repeatable)
  --only         only clear when the named filters are the sole codes
  --passing      only output variants that pass after clearing`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		clearSet := map[string]bool{}
		for _, v := range vcfClearFilterFilters {
			for _, f := range strings.Split(v, ",") {
				clearSet[f] = true
			}
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
		if _, ok := header.InfoDef("CG_CLEARED_FILTER"); !ok {
			header.AddInfo(&vcf.AnnotationDef{
				IsInfo: true, ID: "CG_CLEARED_FILTER", Number: ".", Type: "String",
				Description: "Filters that have been removed from this variant",
			})
		}
		stampVcfProvenance(header, "vcf-clearfilter")

		writer, closeFn, err := openVcfWriter(cmd, vcfClearFilterOutput)
		if err != nil {
			return err
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
			if rec.IsFiltered() {
				clearRecordFilters(rec, clearSet)
			}
			if vcfClearFilterPassing && rec.IsFiltered() {
				continue
			}
			if err := writer.WriteRecord(rec); err != nil {
				return err
			}
		}
		if closeFn != nil {
			return closeFn()
		}
		return writer.Close()
	},
}

// clearRecordFilters removes the targeted FILTER codes from rec and records the
// cleared codes into CG_CLEARED_FILTER, mirroring ngsutilsj's logic.
func clearRecordFilters(rec *vcf.VcfRecord, clearSet map[string]bool) {
	filters := rec.Filters()
	var residual, cleared []string
	if vcfClearFilterOnly {
		only := true
		for _, f := range filters {
			if !clearSet[f] {
				only = false
			}
		}
		if only {
			cleared = append(cleared, filters...)
		} else {
			residual = append(residual, filters...)
		}
	} else {
		for _, f := range filters {
			if clearSet[f] {
				cleared = append(cleared, f)
			} else {
				residual = append(residual, f)
			}
		}
	}
	rec.SetFilters(residual)
	if len(cleared) == 0 {
		return
	}
	val := strings.Join(cleared, ",")
	if existing, ok := rec.Info().Get("CG_CLEARED_FILTER"); ok && existing.String() != "" {
		val = existing.String() + "," + val
	}
	rec.AddInfo("CG_CLEARED_FILTER", val)
}

func init() {
	f := vcfClearFilterCmd.Flags()
	f.StringVarP(&vcfClearFilterOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
	f.StringArrayVar(&vcfClearFilterFilters, "filter", nil, "Filter code to clear (comma-separated, repeatable)")
	f.BoolVar(&vcfClearFilterOnly, "only", false, "Only clear when the named filters are the only codes")
	f.BoolVar(&vcfClearFilterPassing, "passing", false, "Only output passing variants")
}
