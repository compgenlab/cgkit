package seqcmd

import (
	"fmt"

	"github.com/compgenlab/cghts/seqio"
	"github.com/spf13/cobra"
)

// revcompCmd implements the seq-revcomp command: reverse-complement of a sequence.
var revcompCmd = &cobra.Command{
	GroupID:     "seqcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "seq-revcomp seq",
	Short:       "Calculate the reverse-complement of the seq",
	Args:        cobra.ExactArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {

		seq := seqio.NewStringSeq(args[0], "")
		fmt.Fprintln(cmd.OutOrStdout(), seq.FullSeq().RevComp().Seq())

		return nil
	},
}

func init() {
}
