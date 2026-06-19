package bedcmd

import (
	"github.com/compgenlab/hts/bed"
	"github.com/spf13/cobra"
)

var bedToBed6Output string

var bedToBed6Cmd = &cobra.Command{
	GroupID:     "bedcmd",
	Annotations: map[string]string{"since": "v0.3.1"},
	Use:         "bed-tobed6 <input.bed>",
	Short:       "Convert a BED6+ file to a strict BED6 file",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		opts := bed.NewBedWriterOpts().Columns(bed.Columns6)
		return streamBed(cmd, args[0], opts, bedToBed6Output)
	},
}

func init() {
	bedToBed6Cmd.Flags().StringVarP(&bedToBed6Output, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
}
