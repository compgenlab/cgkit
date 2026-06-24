package tabcmd

import (
	"fmt"
	"strings"

	"github.com/compgenlab/hts/htsio/tabix"
	"github.com/spf13/cobra"
)

var tabixIndexCmd = &cobra.Command{
	GroupID:     "tabcmd",
	Annotations: map[string]string{"since": "v0.3.2"},
	Use:         "tabix-index <file.gz>",
	Short:       "Build a tabix (.tbi) index for an existing BGZF-compressed file",
	Long: `Build a tabix .tbi index for a file that is already BGZF-compressed and sorted
by genomic coordinate (the file itself is not modified). The companion index is
written as <file.gz>.tbi.

Use --preset bed, vcf, or gff, or set the columns manually with --seq/--begin/
--end (plus --meta, --skip, --zero-based).`,
	Args: cobra.ExactArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		opts, err := tabixIndexOpts()
		if err != nil {
			return err
		}
		return tabix.NewIndexWriter(opts).WriteIndex(args[0])
	},
}

// tabixIndexOpts builds the tabix configuration from the flags.
func tabixIndexOpts() (*tabix.WriterOpts, error) {
	opts := tabix.NewWriterOpts()
	switch strings.ToLower(tabixIndexPreset) {
	case "bed":
		opts = opts.BED()
	case "vcf":
		opts = opts.VCF()
	case "gff", "gtf":
		opts = opts.GFF()
	case "":
		opts = opts.Columns(tabixIndexColSeq, tabixIndexColBeg, tabixIndexColEnd)
		if tabixIndexZeroBased {
			opts = opts.ZeroBased()
		}
	default:
		return nil, fmt.Errorf("unknown preset %q (use bed, vcf, or gff)", tabixIndexPreset)
	}
	if tabixIndexMeta != "" {
		opts = opts.Meta(tabixIndexMeta[0])
	}
	if tabixIndexSkip > 0 {
		opts = opts.Skip(tabixIndexSkip)
	}
	return opts, nil
}

var (
	tabixIndexPreset    string
	tabixIndexColSeq    int
	tabixIndexColBeg    int
	tabixIndexColEnd    int
	tabixIndexZeroBased bool
	tabixIndexMeta      string
	tabixIndexSkip      int
)

func init() {
	f := tabixIndexCmd.Flags()
	f.StringVarP(&tabixIndexPreset, "preset", "p", "", "Format preset: bed, vcf, gff")
	f.IntVar(&tabixIndexColSeq, "seq", 1, "1-based column for the sequence name")
	f.IntVar(&tabixIndexColBeg, "begin", 2, "1-based column for the start position")
	f.IntVar(&tabixIndexColEnd, "end", 3, "1-based column for the end position (0 for point features)")
	f.BoolVar(&tabixIndexZeroBased, "zero-based", false, "Coordinates are 0-based half-open (like BED)")
	f.StringVar(&tabixIndexMeta, "meta", "", "Comment/header character (e.g. '#')")
	f.IntVar(&tabixIndexSkip, "skip", 0, "Number of header lines to skip")
}
