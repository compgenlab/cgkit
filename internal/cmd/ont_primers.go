package cmd

import (
	"github.com/spf13/cobra"
)

// fastagcCmd implements the initial counting entrypoint.
var ontPrimers = &cobra.Command{
	Use:   "ont-primers <input.fastq>",
	Short: "Find and trim common ONT primers from the start of reads in a FASTQ file",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		return nil
	},
}

func init() {
	rootCmd.AddCommand(ontPrimers)
}
