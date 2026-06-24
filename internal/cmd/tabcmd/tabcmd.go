package tabcmd

import "github.com/spf13/cobra"

func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(tabSortCmd)
	rootCmd.AddCommand(tabixIndexCmd)
	rootCmd.AddGroup(&cobra.Group{ID: "tabcmd", Title: "Tabix"})
}
