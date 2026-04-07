package ontcmd

import "github.com/spf13/cobra"

func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(ontTagsCmd)
	rootCmd.AddCommand(ontUmiClusterCmd)
	rootCmd.AddGroup(&cobra.Group{ID: "ontcmd", Title: "Oxford Nanopore"})
}
