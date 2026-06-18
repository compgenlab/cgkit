package tabcmd

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"github.com/compgenlab/hts/htsio/tabix"
	"github.com/spf13/cobra"
)

var tabSortCmd = &cobra.Command{
	GroupID: "tabcmd",
	Use:     "tab-sort [input] -o <output.bed.gz>",
	Short:   "Sort a tab-delimited file and write as BGZF with optional tabix index",
	Long: `Read a tab-delimited text file (or stdin with '-'), sort by genomic
coordinates, and write a BGZF-compressed output with optional .tbi index.

Preset formats (--preset): bed, vcf, gff. Or specify columns manually
with --seq, --begin, --end.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if tabSortOutput == "" {
			return fmt.Errorf("--output (-o) is required")
		}
		if len(args) < 1 {
			return fmt.Errorf("input file is required (use '-' for stdin)")
		}

		// Build TabixWriter options.
		opts := tabix.NewWriterOpts()
		switch strings.ToLower(tabSortPreset) {
		case "bed":
			opts = opts.BED()
		case "vcf":
			opts = opts.VCF()
		case "gff", "gtf":
			opts = opts.GFF()
		case "":
			opts = opts.Columns(tabSortColSeq, tabSortColBeg, tabSortColEnd)
			if tabSortZeroBased {
				opts = opts.ZeroBased()
			}
		default:
			return fmt.Errorf("unknown preset %q (use bed, vcf, or gff)", tabSortPreset)
		}

		if tabSortMeta != "" {
			opts = opts.Meta(tabSortMeta[0])
		}
		if tabSortSkip > 0 {
			opts = opts.Skip(tabSortSkip)
		}
		if !tabSortNoIndex {
			opts = opts.AutoIndex()
		}

		tw := tabix.NewWriter(tabSortOutput, opts)

		// Open input.
		var scanner *bufio.Scanner
		if args[0] == "-" {
			scanner = bufio.NewScanner(os.Stdin)
		} else {
			f, err := os.Open(args[0])
			if err != nil {
				return fmt.Errorf("opening input: %w", err)
			}
			defer f.Close()
			scanner = bufio.NewScanner(f)
		}
		scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)

		metaCh := byte(0)
		if tabSortMeta != "" {
			metaCh = tabSortMeta[0]
		} else {
			switch strings.ToLower(tabSortPreset) {
			case "vcf", "gff", "gtf":
				metaCh = '#'
			}
		}

		for scanner.Scan() {
			line := scanner.Text()
			if metaCh != 0 && len(line) > 0 && line[0] == metaCh {
				tw.WriteHeader(line)
				continue
			}
			if err := tw.Write(line); err != nil {
				tw.Close()
				return fmt.Errorf("writing line: %w", err)
			}
		}
		if err := scanner.Err(); err != nil {
			tw.Close()
			return fmt.Errorf("reading input: %w", err)
		}

		return tw.Close()
	},
}

var (
	tabSortOutput    string
	tabSortPreset    string
	tabSortColSeq    int
	tabSortColBeg    int
	tabSortColEnd    int
	tabSortZeroBased bool
	tabSortMeta      string
	tabSortSkip      int
	tabSortNoIndex   bool
)

func init() {
	tabSortCmd.Flags().StringVarP(&tabSortOutput, "output", "o", "", "Output BGZF file path (required)")
	tabSortCmd.Flags().StringVarP(&tabSortPreset, "preset", "p", "", "Format preset: bed, vcf, gff")
	tabSortCmd.Flags().IntVar(&tabSortColSeq, "seq", 1, "1-based column for sequence name")
	tabSortCmd.Flags().IntVar(&tabSortColBeg, "begin", 2, "1-based column for start position")
	tabSortCmd.Flags().IntVar(&tabSortColEnd, "end", 3, "1-based column for end position (0 for point features)")
	tabSortCmd.Flags().BoolVar(&tabSortZeroBased, "zero-based", false, "Coordinates are 0-based half-open (like BED)")
	tabSortCmd.Flags().StringVar(&tabSortMeta, "meta", "", "Comment/header character (e.g. '#')")
	tabSortCmd.Flags().IntVar(&tabSortSkip, "skip", 0, "Number of header lines to skip")
	tabSortCmd.Flags().BoolVar(&tabSortNoIndex, "no-index", false, "Disable automatic .tbi index generation")
}
