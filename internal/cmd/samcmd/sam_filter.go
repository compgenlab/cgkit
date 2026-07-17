package samcmd

import (
	"fmt"
	"strings"

	"github.com/compgenlab/cghts/htsio"
	"github.com/compgenlab/cghts/htsio/bam"
	"github.com/compgenlab/cghts/htsio/cram"
	"github.com/compgenlab/cghts/htsio/sam"
	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/spf13/cobra"
)

var samFilterCmd = &cobra.Command{
	GroupID:     "samcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "sam-filter <input.bam> [output]",
	Short:       "Filter SAM/BAM/CRAM reads and write to a new file",
	Long: "Filter reads from a SAM/BAM/CRAM file and write passing reads to a new file.\n\n" +
		"The output format defaults to \"auto\": when an output file is given, the format\n" +
		"is chosen from its extension (.bam → BAM, .cram → CRAM, otherwise SAM text).\n" +
		"With no output file, SAM text is written to stdout. Use --format to force a\n" +
		"specific format regardless of the output filename.",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) < 1 {
			cmd.Help()
			return nil
		}

		opts, err := samFilterReaderFlags.buildReaderOpts()
		if err != nil {
			return err
		}

		inputFile := args[0]
		outputFile := ""
		if len(args) >= 2 {
			outputFile = args[1]
		}

		format, err := resolveSamFilterFormat(outputFile)
		if err != nil {
			return err
		}

		reader, err := htsio.NewSamReader(inputFile, opts)
		if err != nil {
			return err
		}
		defer reader.Close()

		header, err := reader.Header()
		if err != nil {
			return fmt.Errorf("reading header: %w", err)
		}
		header.AddPGLine("sam-filter", "cgkit", buildinfo.String())

		// Create the writer for the resolved output format.
		var writer htsio.SamWriter
		switch format {
		case "bam":
			w, err := bam.NewWriter(outputFile, header)
			if err != nil {
				return err
			}
			writer = w
		case "cram":
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
		default: // "sam"
			// Honor the output filename for SAM too; stdout when unspecified.
			samOut := outputFile
			if samOut == "" {
				samOut = "-"
			}
			w, err := sam.NewWriter(samOut, header)
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
	samFilterFormat      string
	samFilterReaderFlags samReaderFlags
)

// resolveSamFilterFormat determines the output format ("sam", "bam", or "cram")
// from the --format flag and the output filename.
//
// The default format is "auto": it is inferred from the output file extension
// (.bam → BAM, .cram → CRAM, anything else → SAM text), falling back to SAM
// text on stdout when no output file is given.
func resolveSamFilterFormat(outputFile string) (string, error) {
	switch strings.ToLower(samFilterFormat) {
	case "sam", "bam", "cram":
		return strings.ToLower(samFilterFormat), nil
	case "", "auto":
		if outputFile == "" || outputFile == "-" {
			return "sam", nil
		}
		lower := strings.ToLower(outputFile)
		switch {
		case strings.HasSuffix(lower, ".bam"):
			return "bam", nil
		case strings.HasSuffix(lower, ".cram"):
			return "cram", nil
		default:
			return "sam", nil
		}
	default:
		return "", fmt.Errorf("unknown output format %q (expected auto, sam, bam, or cram)", samFilterFormat)
	}
}

func init() {
	samFilterReaderFlags.register(samFilterCmd)

	samFilterCmd.Flags().StringVarP(&samFilterFormat, "format", "O", "auto", "Output format: auto, sam, bam, or cram")
}
