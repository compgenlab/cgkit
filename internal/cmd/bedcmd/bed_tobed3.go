package bedcmd

import (
	"github.com/compgenlab/cghts/bed"
	"github.com/spf13/cobra"
)

var bedToBed3Output string

var bedToBed3Cmd = &cobra.Command{
	GroupID:     "bedcmd",
	Annotations: map[string]string{"since": "v0.3.1"},
	Use:         "bed-tobed3 <input.bed>",
	Short:       "Convert a BED3+ file to a strict BED3 file",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		opts := bed.NewBedWriterOpts().Columns(bed.Columns3)
		return streamBed(cmd, args[0], opts, bedToBed3Output)
	},
}

func init() {
	bedToBed3Cmd.Flags().StringVarP(&bedToBed3Output, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
}
