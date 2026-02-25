package ontcmd

import "github.com/spf13/cobra"

func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(ontPrimersCmd)
	rootCmd.AddGroup(&cobra.Group{ID: "ontcmd", Title: "Oxford Nanopore"})

}
