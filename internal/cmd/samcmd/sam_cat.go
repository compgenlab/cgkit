package samcmd

import (
	"fmt"
	"io"

	"github.com/compgen-io/cgltk/htsio"
	"github.com/spf13/cobra"
)

var samCatCmd = &cobra.Command{
	GroupID: "samcmd",
	Use:     "sam-cat <input>",
	Short:   "Read a SAM/BAM/CRAM file and write SAM text to stdout",
	Hidden:  true,
	Args:    cobra.ExactArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		reader, err := htsio.NewSamReader(args[0], nil)
		if err != nil {
			return fmt.Errorf("open %s: %w", args[0], err)
		}
		defer reader.Close()

		header, err := reader.Header()
		if err != nil {
			return fmt.Errorf("read header: %w", err)
		}

		writer := htsio.NewStdoutSamWriter(header)
		defer writer.Close()

		for {
			rec, err := reader.Next()
			if err == io.EOF {
				break
			}
			if err != nil {
				return fmt.Errorf("read record: %w", err)
			}
			if err := writer.Write(rec); err != nil {
				return fmt.Errorf("write record: %w", err)
			}
		}
		return nil
	},
}
