package fastqcmd

import (
	"fmt"
	"io"

	"github.com/compgenlab/hts/seqio"
	"github.com/spf13/cobra"
)

var fastqTagCmd = &cobra.Command{
	GroupID:     "fastaqcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "fastq-tag <tag> <input.fastq>",
	Short:       "Add a tag to the comment field of FASTQ records",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) < 2 {
			cmd.Help()
			return nil
		}
		tag := args[0]
		reader, err := seqio.NewFastqFile(args[1])
		if err != nil {
			return err
		}
		defer reader.Close()

		for rec, err := reader.NextFastqSeq(); ; rec, err = reader.NextFastqSeq() {
			if err != nil {
				if err == io.EOF {
					break
				}
				return err
			}
			sq := rec.FullSeq()
			comment := rec.Comment()
			if comment != "" {
				comment = comment + "\t" + tag
			} else {
				comment = tag
			}
			fmt.Fprintf(cmd.OutOrStdout(), "@%s %s\n%s\n+\n%s\n", rec.Name(), comment, sq.Seq(), sq.Qual())
		}
		return nil
	},
}
