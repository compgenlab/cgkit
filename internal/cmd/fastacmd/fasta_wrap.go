package fastacmd

import (
	"fmt"
	"io"

	"github.com/compgen-io/cgltk/seqio"
	"github.com/spf13/cobra"
)

// fastagcCmd implements the initial counting entrypoint.
var fastaWrapCmd = &cobra.Command{
	GroupID: "fastacmd",
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
			if wrapWidth < 1 {
				// if wrapWidth is negative, just print the whole sequence at once
				for chunk := range rec.Chunks(1024) {
					fmt.Print(chunk.Seq())
				}
				fmt.Print('\n')
				continue
			} else {
				seq := ""
				for chunk := range rec.Chunks(wrapWidth) {
					seq += chunk.Seq()
					for len(seq) >= wrapWidth {
						fmt.Printf("%s\n", seq[:wrapWidth])
						seq = seq[wrapWidth:]
					}
				}
				for len(seq) >= wrapWidth {
					fmt.Printf("%s\n", seq[:wrapWidth])
					seq = seq[wrapWidth:]
				}
				if len(seq) > 0 {
					fmt.Printf("%s\n", seq)
				}
			}
		}
		return nil
	},
}

var wrapWidth int

func init() {
	fastaWrapCmd.Flags().IntVarP(&wrapWidth, "width", "w", 70, "Line width to wrap sequences to (-1 for no wrapping)")
}
