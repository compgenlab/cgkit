package seqcmd

import "github.com/spf13/cobra"

func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(swalignCmd)
	rootCmd.AddCommand(revcompCmd)
	rootCmd.AddCommand(msaCmd)

	rootCmd.AddGroup(&cobra.Group{ID: "seqcmd", Title: "Sequence"})
}
