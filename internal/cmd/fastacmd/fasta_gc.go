package fastacmd

import (
	"io"

	seqanalysis "github.com/compgenlab/cghts/analysis/seq"
	"github.com/compgenlab/cghts/seqio"
	"github.com/spf13/cobra"
)

// fastaGCCmd implements the fasta-gc command: per-sequence GC content.
var fastaGCCmd = &cobra.Command{
	GroupID:     "fastaqcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "fasta-gc <input.fasta>",
	Short:       "Return the GC content of sequences in a FASTA file",
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
			pct := seqanalysis.CalcGC(rec)
			cmd.Printf("%s\t%.4f\n", rec.Name(), pct)
		}
		return nil
	},
}
