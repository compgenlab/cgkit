package fastqcmd

import "github.com/spf13/cobra"

func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(fastqGCCmd)
	rootCmd.AddCommand(fastqTagCmd)
	rootCmd.AddGroup(&cobra.Group{ID: "fastqcmd", Title: "FASTQ"})

}
