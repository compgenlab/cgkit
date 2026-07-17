package fastqcmd

import (
	"io"

	seqanalysis "github.com/compgenlab/cghts/analysis/seq"
	"github.com/compgenlab/cghts/seqio"
	"github.com/spf13/cobra"
)

// fastqGCCmd implements the fastq-gc command: per-sequence GC content.
var fastqGCCmd = &cobra.Command{
	GroupID:     "fastaqcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "fastq-gc <input.fastq>",
	Short:       "Return the GC content of sequences in a FASTQ file",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		reader, err := seqio.NewFastqFile(args[0])
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
