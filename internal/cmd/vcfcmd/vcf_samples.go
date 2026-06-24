package vcfcmd

import (
	"fmt"

	"github.com/spf13/cobra"
)

var vcfSamplesCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-samples <input.vcf>",
	Short:       "Output the sample names in a VCF file",
	Long:        "Output the sample names in a VCF file, one per line. Use - to read from stdin.",
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
		out := cmd.OutOrStdout()
		for _, sample := range header.Samples() {
			fmt.Fprintln(out, sample)
		}
		return nil
	},
}
