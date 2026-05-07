package seqcmd

import (
	"fmt"

	"github.com/compgen-io/cgkit/seqio"
	"github.com/spf13/cobra"
)

// fastagcCmd implements the initial counting entrypoint.
var revcompCmd = &cobra.Command{
	GroupID: "seqcmd",
	Use:     "seq-revcomp seq",
	Short:   "Calculate the reverse-compliment of the seq",
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
