package cmd

import (
	"fmt"
	"os"

	"github.com/compgen-io/cgltk/internal/cmd/fastacmd"
	"github.com/compgen-io/cgltk/internal/cmd/fastqcmd"
	"github.com/compgen-io/cgltk/internal/cmd/ontcmd"
	"github.com/compgen-io/cgltk/internal/cmd/samcmd"
	"github.com/compgen-io/cgltk/internal/cmd/seqcmd"
	"github.com/spf13/cobra"
)

// Set via -ldflags at build time.
var (
	Version = "dev"
	GitHash = ""
)

func versionString() string {
	if GitHash != "" {
		return fmt.Sprintf("%s (%s)", Version, GitHash)
	}
	return Version
}

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:     "cgltk",
	Short:   "Toolkit for computational genomics research",
	Version: "", // set in init()
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
	rootCmd.Version = versionString()
	// SilenceUsage prevents Cobra from printing usage on errors after argument parsing.
	rootCmd.SilenceUsage = true
	rootCmd.CompletionOptions.HiddenDefaultCmd = true

	// Add version footer to all help output.
	defaultHelp := rootCmd.HelpTemplate()
	rootCmd.SetHelpTemplate(defaultHelp + "\ncgltk " + versionString() + "\n")

	ontcmd.InitCmd(rootCmd)
	fastacmd.InitCmd(rootCmd)
	fastqcmd.InitCmd(rootCmd)
	samcmd.InitCmd(rootCmd)
	seqcmd.InitCmd(rootCmd)
}
