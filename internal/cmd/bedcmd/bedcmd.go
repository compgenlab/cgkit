package bedcmd

import "github.com/spf13/cobra"

func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(bedCleanCmd)
	rootCmd.AddCommand(bedToBed3Cmd)
	rootCmd.AddCommand(bedToBed6Cmd)
	rootCmd.AddCommand(bedResizeCmd)
	rootCmd.AddCommand(bedStatsCmd)
	rootCmd.AddCommand(bedToFastaCmd)
	rootCmd.AddCommand(bedSetCmd)
	rootCmd.AddGroup(&cobra.Group{ID: "bedcmd", Title: "BED"})
}
