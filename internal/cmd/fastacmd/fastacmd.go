package fastacmd

import "github.com/spf13/cobra"

func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(fastaCatCmd)
	rootCmd.AddCommand(fastaWrapCmd)
	rootCmd.AddCommand(fastaGCCmd)
	rootCmd.AddGroup(&cobra.Group{ID: "fastacmd", Title: "FASTA"})

}
