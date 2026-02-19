package cmd

import (
	"io"

	"github.com/compgen-io/cgltk/seqio"
	"github.com/compgen-io/cgltk/sequtils"
	"github.com/spf13/cobra"
)

// fastagcCmd implements the initial counting entrypoint.
var fastaGCCmd = &cobra.Command{
	Use:   "fasta-gc <input.fasta>",
	Short: "Return the GC content of sequences in a FASTA file",
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
			if err != nil && err != io.EOF {
				return err
			}
			if rec == nil {
				break
			}
			pct := sequtils.CalcGC(rec)
			cmd.Printf("%s\t%.4f\n", rec.Name(), pct)
		}
		return nil
	},
}

func init() {
	rootCmd.AddCommand(fastaGCCmd)
}
