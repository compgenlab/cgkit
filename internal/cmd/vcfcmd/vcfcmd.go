package vcfcmd

import "github.com/spf13/cobra"

// InitCmd registers the VCF command group with the root command.
func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(vcfSamplesCmd)
	rootCmd.AddCommand(vcfToBedCmd)
	rootCmd.AddCommand(vcfExportCmd)
	rootCmd.AddCommand(vcfReorderCmd)
	rootCmd.AddCommand(vcfAnnotateCmd)
	rootCmd.AddCommand(vcfFilterCmd)
	rootCmd.AddCommand(vcfStatsCmd)
	rootCmd.AddCommand(vcfTsTvCmd)
	rootCmd.AddCommand(vcfClearFilterCmd)
	rootCmd.AddCommand(vcfRenameCmd)
	rootCmd.AddCommand(vcfChrFixCmd)
	rootCmd.AddCommand(vcfRemoveFlagsCmd)
	rootCmd.AddCommand(vcfHeaderInfoCmd)
	rootCmd.AddCommand(vcfCheckCmd)
	rootCmd.AddCommand(vcfSplitCmd)
	rootCmd.AddCommand(vcfSampleExportCmd)
	rootCmd.AddCommand(vcfToCountCmd)
	rootCmd.AddCommand(vcfStripCmd)
	rootCmd.AddCommand(vcfConcatCmd)
	rootCmd.AddCommand(vcfMergeCmd)
	rootCmd.AddGroup(&cobra.Group{ID: "vcfcmd", Title: "VCF"})
}
