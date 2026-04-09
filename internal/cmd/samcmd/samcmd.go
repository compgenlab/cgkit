package samcmd

import "github.com/spf13/cobra"

func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(samExportCmd)
	rootCmd.AddCommand(samToSeqCmd)
	rootCmd.AddGroup(&cobra.Group{ID: "samcmd", Title: "SAM/BAM/CRAM"})
}
