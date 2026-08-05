package bedcmd

import (
	"github.com/compgenlab/cghts/bed"
	"github.com/spf13/cobra"
)

var bedCleanOutput string

var bedCleanCmd = &cobra.Command{
	GroupID:     "bedcmd",
	Annotations: map[string]string{"since": "v0.3.1"},
	Use:         "bed-clean <input.bed>",
	Short:       "Clean BED score entries to be integers (expands records to BED6+)",
	Args:        cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		opts := bed.NewBedWriterOpts().ForceScoreInt(true)
		return streamBed(cmd, args[0], opts, bedCleanOutput)
	},
}

func init() {
	bedCleanCmd.Flags().StringVarP(&bedCleanOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
}
