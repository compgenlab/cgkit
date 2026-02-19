package cmd

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "cgltk",
	Short: "Toolkit for computational genomics research",
	Long: `Utility toolkit for computational genomics research, 
with a collection of commands for sequence analysis, 
NGS data-wrangling, and more.`,
}

// Execute runs the root command and exits with a non-zero status on error.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func init() {
	// SilenceUsage prevents Cobra from printing usage on errors after argument parsing.
	rootCmd.SilenceUsage = true
	rootCmd.CompletionOptions.HiddenDefaultCmd = true
}
