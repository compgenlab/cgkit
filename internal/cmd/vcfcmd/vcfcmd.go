package vcfcmd

import "github.com/spf13/cobra"

// InitCmd registers the VCF command group with the root command.
func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(vcfSamplesCmd)
	rootCmd.AddCommand(vcfToBedCmd)
	rootCmd.AddCommand(vcfExportCmd)
	rootCmd.AddCommand(vcfReorderCmd)
	rootCmd.AddCommand(vcfStatsCmd)
	rootCmd.AddCommand(vcfTsTvCmd)
	rootCmd.AddGroup(&cobra.Group{ID: "vcfcmd", Title: "VCF"})
}
