package cmd

import (
	"fmt"
	"os"

	"github.com/compgenlab/cgio/internal/buildinfo"
	"github.com/compgenlab/cgio/internal/cmd/bedcmd"
	"github.com/compgenlab/cgio/internal/cmd/fastacmd"
	"github.com/compgenlab/cgio/internal/cmd/fastqcmd"
	"github.com/compgenlab/cgio/internal/cmd/ontcmd"
	"github.com/compgenlab/cgio/internal/cmd/samcmd"
	"github.com/compgenlab/cgio/internal/cmd/seqcmd"
	"github.com/compgenlab/cgio/internal/cmd/tabcmd"
	"github.com/compgenlab/cgio/internal/cmd/vcfcmd"
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
	Use:     "cgio",
	Short:   "Toolkit for computational genomics research",
	Version: "", // set in init()
	Long: `Utility toolkit for computational genomics research,
with a collection of commands for sequence analysis,
NGS data-wrangling, and more.`,
	// Silence usage only once we reach a command's RunE. Cobra validates args
	// and flags *before* PersistentPreRun, so parsing errors (e.g. a missing
	// argument) still print usage, while errors during execution do not.
	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		cmd.SilenceUsage = true
	},
}

// Execute runs the root command and exits with a non-zero status on error.
// Cobra already prints the error (and usage, for argument/flag errors), so we
// only need to set the exit status here.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		os.Exit(1)
	}
}

func init() {
	// Share the ldflags-injected version with the leaf buildinfo package so
	// command groups can stamp provenance without importing this package.
	buildinfo.Version = Version
	buildinfo.GitHash = GitHash

	rootCmd.Version = versionString()
	rootCmd.CompletionOptions.HiddenDefaultCmd = true

	// Help template: cobra's default command help, plus a "Since:" line sourced
	// from each command's "since" annotation, plus a version footer.
	const helpBody = `{{with (or .Long .Short)}}{{. | trimTrailingWhitespaces}}

{{end}}{{if .Annotations}}{{with (index .Annotations "since")}}Since: {{.}}

{{end}}{{end}}{{if or .Runnable .HasSubCommands}}{{.UsageString}}{{end}}`
	rootCmd.SetHelpTemplate(helpBody + "\ncgio " + versionString() +
		" https://compgenlab.org/cgio\n")

	ontcmd.InitCmd(rootCmd)
	bedcmd.InitCmd(rootCmd)
	fastacmd.InitCmd(rootCmd)
	fastqcmd.InitCmd(rootCmd)
	samcmd.InitCmd(rootCmd)
	seqcmd.InitCmd(rootCmd)
	tabcmd.InitCmd(rootCmd)
	vcfcmd.InitCmd(rootCmd)
}
