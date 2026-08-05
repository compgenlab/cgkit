package fastacmd

import "github.com/spf13/cobra"

func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(fastaCatCmd)
	rootCmd.AddCommand(fastaWrapCmd)
	rootCmd.AddCommand(fastaGCCmd)
	// Guarded because fastqcmd registers the same shared group, and
	// whichever package is initialized first should win rather than the
	// second one adding a duplicate.
	if !rootCmd.ContainsGroup("fastaqcmd") {
		rootCmd.AddGroup(&cobra.Group{ID: "fastaqcmd", Title: "FASTA/Q"})
	}
}
