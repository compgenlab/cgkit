package fastacmd

import (
	"io"
	"os"

	"github.com/compgen-io/cgltk/seqio"
	"github.com/spf13/cobra"
)

var fastaWrapCmd = &cobra.Command{
	GroupID: "fastaqcmd",
	Use:     "fasta-wrap <input.fasta>",
	Short:   "Reformat the sequences in a FASTA file to a specified line width",
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

		wrap := wrapWidth
		if wrap < 1 {
			wrap = 0
		}
		writer := seqio.NewFastaWriter(os.Stdout, seqio.NewFastaWriterOpts().Wrap(wrap))
		defer writer.Close()

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

			if err := writer.WriteSeq(rec); err != nil {
				return err
			}
		}
		return nil
	},
}

var wrapWidth int

func init() {
	fastaWrapCmd.Flags().IntVarP(&wrapWidth, "width", "w", 70, "Line width to wrap sequences to (-1 for no wrapping)")
}
