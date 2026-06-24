package samcmd

import (
	"fmt"

	"github.com/compgenlab/cgio/internal/buildinfo"
	"github.com/compgenlab/hts/htsio"
	"github.com/compgenlab/hts/htsio/bam"
	"github.com/compgenlab/hts/htsio/cram"
	"github.com/compgenlab/hts/htsio/sam"
	"github.com/spf13/cobra"
)

var samFilterCmd = &cobra.Command{
	GroupID:     "samcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "sam-filter <input.bam> <output.bam>",
	Short:       "Filter SAM/BAM/CRAM reads and write to a new file",
	Long:        "Filter reads from a SAM/BAM/CRAM file and write passing reads to a new SAM/BAM/CRAM file.",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) < 2 {
			cmd.Help()
			return nil
		}

		opts, err := samFilterReaderFlags.buildReaderOpts()
		if err != nil {
			return err
		}

		inputFile := args[0]
		outputFile := args[1]

		reader, err := htsio.NewSamReader(inputFile, opts)
		if err != nil {
			return err
		}
		defer reader.Close()

		header, err := reader.Header()
		if err != nil {
			return fmt.Errorf("reading header: %w", err)
		}
		header.AddPGLine("sam-filter", "cgio", buildinfo.String())

		// Determine output format.
		var writer htsio.SamWriter
		if samFilterBAM {
			w, err := bam.NewWriter(outputFile, header)
			if err != nil {
				return err
			}
			writer = w
		} else if samFilterCRAM {
			// Reuse --cram-ref as the output reference when provided;
			// otherwise write a reference-free CRAM (literal bases).
			cramOpts := cram.NewWriterOpts()
			if ref := samFilterReaderFlags.cramRef; ref != "" {
				cramOpts.Reference(ref)
			}
			w, err := cram.NewWriter(outputFile, header, cramOpts)
			if err != nil {
				return err
			}
			writer = w
		} else {
			w, err := sam.NewWriter("-", header)
			if err != nil {
				return err
			}
			writer = w
		}
		defer writer.Close()

		// If a region is specified, use Query() to iterate the subset.
		if regionStr := samFilterReaderFlags.queryRegion(); regionStr != "" {
			ref, start, end, err := htsio.ParseRegion(regionStr)
			if err != nil {
				return err
			}
			if end < 0 {
				end = 1<<30 - 1
			}
			records, err := reader.Query(ref, start, end)
			if err != nil {
				return fmt.Errorf("query %q: %w", regionStr, err)
			}
			for rec, err := range records {
				if err != nil {
					return err
				}
				if err := writer.Write(rec); err != nil {
					return fmt.Errorf("write record: %w", err)
				}
			}
			return nil
		}

		// No region — stream all records.
		for rec, err := range reader.Records() {
			if err != nil {
				return err
			}
			if err := writer.Write(rec); err != nil {
				return fmt.Errorf("write record: %w", err)
			}
		}

		return nil
	},
}

var (
	samFilterBAM         bool
	samFilterCRAM        bool
	samFilterReaderFlags samReaderFlags
)

func init() {
	samFilterReaderFlags.register(samFilterCmd)

	samFilterCmd.Flags().BoolVar(&samFilterBAM, "bam", false, "Output in BAM format")
	samFilterCmd.Flags().BoolVar(&samFilterCRAM, "cram", false, "Output in CRAM format")
}
