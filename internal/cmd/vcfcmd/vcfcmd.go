package vcfcmd

import "github.com/spf13/cobra"

// InitCmd registers the VCF command group with the root command.
func InitCmd(rootCmd *cobra.Command) {
	vcfVarstoreCmd.AddCommand(vcfVarstoreCreateCmd)
	vcfVarstoreCmd.AddCommand(vcfVarstoreAttachCmd)
	rootCmd.AddCommand(vcfVarstoreCmd)
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
	rootCmd.AddCommand(vcfGtCountCmd)
	rootCmd.AddCommand(vcfVarSummaryCmd)
	rootCmd.AddCommand(vcfStripCmd)
	rootCmd.AddCommand(vcfConcatCmd)
	rootCmd.AddCommand(vcfMergeCmd)
	rootCmd.AddCommand(vcfToBedpeCmd)
	rootCmd.AddCommand(vcfSvToFastaCmd)
	rootCmd.AddCommand(vcfToVarstoreCmd)
	rootCmd.AddCommand(vcfRefBlocksCmd)
	rootCmd.AddCommand(vcfVarQueryCmd)
	rootCmd.AddGroup(&cobra.Group{ID: "vcfcmd", Title: "VCF"})
	// The genotype store has its own section: it is a product with its own
	// pipeline (convert, refblocks, attach, query), not a VCF utility.
	rootCmd.AddGroup(&cobra.Group{ID: "varstore", Title: "Variant store"})
}
