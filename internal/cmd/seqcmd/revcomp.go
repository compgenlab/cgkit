package seqcmd

import (
	"fmt"

	"github.com/compgenlab/hts/seqio"
	"github.com/spf13/cobra"
)

// revcompCmd implements the seq-revcomp command: reverse-complement of a sequence.
var revcompCmd = &cobra.Command{
	GroupID:     "seqcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "seq-revcomp seq",
	Short:       "Calculate the reverse-complement of the seq",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) != 1 {
			cmd.Help()
			return nil
		}

		seq := seqio.NewStringSeq(args[0], "")
		fmt.Println(seq.FullSeq().RevComp().Seq())

		return nil
	},
}

func init() {
}
