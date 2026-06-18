package samcmd

import (
	"fmt"

	"github.com/compgenlab/hts/htsio"
	_ "github.com/compgenlab/hts/htsio/bam"
	_ "github.com/compgenlab/hts/htsio/cram"
	"github.com/compgenlab/hts/htsio/sam"
	"github.com/spf13/cobra"
)

var samCatCramRef string

func init() {
	samCatCmd.Flags().StringVar(&samCatCramRef, "cram-ref", "", "Reference FASTA for CRAM files")
}

var samCatCmd = &cobra.Command{
	GroupID: "samcmd",
	Use:     "sam-cat <input>",
	Short:   "Read a SAM/BAM/CRAM file and write SAM text to stdout",
	Hidden:  true,
	Args:    cobra.ExactArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		opts := htsio.NewSamReaderOpts()
		if samCatCramRef != "" {
			opts.RefPath(samCatCramRef)
		}
		reader, err := htsio.NewSamReader(args[0], opts)
		if err != nil {
			return fmt.Errorf("open %s: %w", args[0], err)
		}
		defer reader.Close()

		header, err := reader.Header()
		if err != nil {
			return fmt.Errorf("read header: %w", err)
		}

		writer, err := sam.NewWriter("-", header)
		if err != nil {
			return fmt.Errorf("create writer: %w", err)
		}
		defer writer.Close()

		for rec, err := range reader.Records() {
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
