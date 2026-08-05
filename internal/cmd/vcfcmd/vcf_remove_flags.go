package vcfcmd

import (
	"fmt"
	"strings"

	"github.com/compgenlab/cghts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfRemoveFlagsOutput string
	vcfRemoveFlagsKey    string
	vcfRemoveFlagsAlways bool
)

var vcfRemoveFlagsCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-remove-flags <input.vcf>",
	Short:       "Replace all INFO flags with a comma-separated list",
	Long: `Replace every Flag-typed INFO field with a single key holding the set flags as
a comma-separated list (FOO;BAR => FLAGS=FOO,BAR).

  --key NAME   name for the new INFO key (default FLAGS)
  --always     always emit the key, using "." when no flags are set`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		flagSet := map[string]bool{}
		return runVcfStream(cmd, vcfStream{
			name: "vcf-remove-flags",
			in:   args[0],
			out:  vcfRemoveFlagsOutput,
			header: func(h *vcf.VcfHeader) error {
				var flagIDs []string
				for _, id := range h.InfoIDs() {
					if d, ok := h.InfoDef(id); ok && strings.EqualFold(d.Type, "Flag") {
						flagIDs = append(flagIDs, id)
						flagSet[id] = true
					}
				}
				if len(flagIDs) == 0 {
					return fmt.Errorf("no INFO flags defined in VCF")
				}
				for _, id := range flagIDs {
					h.RemoveInfo(id)
				}
				h.AddInfo(&vcf.AnnotationDef{
					IsInfo: true, ID: vcfRemoveFlagsKey, Number: ".", Type: "String",
					Description: "INFO Flag values as CSV (" + strings.Join(flagIDs, ",") + ")",
				})
				return nil
			},
			record: func(rec *vcf.VcfRecord) (bool, error) {
				var setFlags []string
				for _, id := range rec.Info().Keys() {
					if flagSet[id] {
						setFlags = append(setFlags, id)
					}
				}
				if len(setFlags) > 0 {
					for _, id := range setFlags {
						rec.Info().Remove(id)
					}
					rec.AddInfo(vcfRemoveFlagsKey, strings.Join(setFlags, ","))
				} else if vcfRemoveFlagsAlways {
					rec.AddInfo(vcfRemoveFlagsKey, ".")
				}
				return true, nil
			},
		})
	},
}

func init() {
	f := vcfRemoveFlagsCmd.Flags()
	addVcfOutputFlags(vcfRemoveFlagsCmd, &vcfRemoveFlagsOutput)
	f.StringVar(&vcfRemoveFlagsKey, "key", "FLAGS", "Name for the new INFO key")
	f.BoolVar(&vcfRemoveFlagsAlways, "always", false, "Always include the key (use '.' when no flags set)")
}
