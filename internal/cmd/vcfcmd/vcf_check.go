package vcfcmd

import (
	"fmt"
	"io"

	"github.com/spf13/cobra"
)

var vcfCheckQuiet bool

var vcfCheckCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-check <input.vcf>",
	Short:       "Validate a VCF file",
	Long: `Validate a VCF file by parsing every record. The command exits with an error
on the first record that fails to parse, and succeeds silently otherwise.`,
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

		if _, err := reader.Header(); err != nil {
			return err
		}
		var count int64
		for {
			_, err := reader.NextRecord()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			count++
		}
		if !vcfCheckQuiet {
			fmt.Fprintf(cmd.ErrOrStderr(), "OK: %d variants\n", count)
		}
		return nil
	},
}

func init() {
	vcfCheckCmd.Flags().BoolVarP(&vcfCheckQuiet, "quiet", "q", false, "Quiet output (no summary)")
}
