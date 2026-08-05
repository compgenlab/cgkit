package vcfcmd

import (
	"strings"

	"github.com/compgenlab/cghts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfClearFilterOutput  string
	vcfClearFilterFilters []string
	vcfClearFilterOnly    bool
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
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		clearSet := map[string]bool{}
		for _, v := range vcfClearFilterFilters {
			for _, f := range strings.Split(v, ",") {
				clearSet[f] = true
			}
		}

		return runVcfStream(cmd, vcfStream{
			name: "vcf-clearfilter",
			in:   args[0],
			out:  vcfClearFilterOutput,
			header: func(h *vcf.VcfHeader) error {
				if _, ok := h.InfoDef("CG_CLEARED_FILTER"); !ok {
					h.AddInfo(&vcf.AnnotationDef{
						IsInfo: true, ID: "CG_CLEARED_FILTER", Number: ".", Type: "String",
						Description: "Filters that have been removed from this variant",
					})
				}
				return nil
			},
			record: func(rec *vcf.VcfRecord) (bool, error) {
				if rec.IsFiltered() {
					clearRecordFilters(rec, clearSet)
				}
				return !(vcfPassing && rec.IsFiltered()), nil
			},
		})
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
	addVcfOutputFlags(vcfClearFilterCmd, &vcfClearFilterOutput)
	f.StringArrayVar(&vcfClearFilterFilters, "filter", nil, "Filter code to clear (comma-separated, repeatable)")
	f.BoolVar(&vcfClearFilterOnly, "only", false, "Only clear when the named filters are the only codes")
	addPassingFlag(vcfClearFilterCmd, "Only output passing variants")
}
