package fastqcmd

import "github.com/spf13/cobra"

func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(fastqGCCmd)
	rootCmd.AddCommand(fastqTagCmd)
	// Group "fastaqcmd" is registered by fastacmd.InitCmd.
}
