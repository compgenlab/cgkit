package fastacmd

import (
	"fmt"
	"io"

	"github.com/compgenlab/hts/seqio"
	"github.com/spf13/cobra"
)

// fastaCatCmd implements the fasta-cat command: read and re-emit FASTA records unwrapped.
var fastaCatCmd = &cobra.Command{
	GroupID:     "fastaqcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "fasta-cat <input.fasta>",
	Short:       "Write the sequences in a FASTA file without any wrapping",
	Hidden:      true,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		reader, err := seqio.NewFastaFile(args[0])
		if err != nil {
			return err
		}
		defer reader.Close()

		for rec, err := reader.NextSeq(); ; rec, err = reader.NextSeq() {
			if err != nil {
				if err != io.EOF {
					return err
				}
				break
			}
			if rec == nil {
				break
			}

			fmt.Printf(">%s", rec.Name())

			if rec.Comment() != "" {
				fmt.Printf(" %s\n", rec.Comment())
			} else {
				fmt.Printf("\n")
			}
			fmt.Printf("%s\n", rec.FullSeq().Seq())
		}
		return nil
	},
}
