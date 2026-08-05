package fastqcmd

import "github.com/spf13/cobra"

func InitCmd(rootCmd *cobra.Command) {
	// FASTA and FASTQ share one help group, and either package may be
	// initialized first. This used to be left to fastacmd, on the strength of
	// root.go calling it earlier -- so fastqcmd.InitCmd on its own registered
	// commands against a group id that did not exist, which cobra panics on.
	// That made the first fastqcmd-only test impossible to write without
	// knowing the trick, and fastq_tag_test.go duly hand-rolled the group.
	if !rootCmd.ContainsGroup("fastaqcmd") {
		rootCmd.AddGroup(&cobra.Group{ID: "fastaqcmd", Title: "FASTA/Q"})
	}
	rootCmd.AddCommand(fastqGCCmd)
	rootCmd.AddCommand(fastqTagCmd)
}
