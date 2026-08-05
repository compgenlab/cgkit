package tabcmd

import (
	"fmt"
	"strings"

	"github.com/compgenlab/cghts/htsio/tabix"
	"github.com/spf13/cobra"
)

func InitCmd(rootCmd *cobra.Command) {
	rootCmd.AddCommand(tabSortCmd)
	rootCmd.AddCommand(tabixIndexCmd)
	rootCmd.AddGroup(&cobra.Group{ID: "tabcmd", Title: "Tabix"})
}

// tabixSpec is the index configuration tab-sort and tabix-index both build.
//
// The two register the same seven flags against separate variables and then ran
// the same preset switch over them, twice, character for character. The flags
// stay separate -- they belong to different commands -- but the interpretation
// is one function, so a preset added to one cannot go missing from the other.
type tabixSpec struct {
	preset                 string
	colSeq, colBeg, colEnd int
	zeroBased              bool
	meta                   string
	skip                   int
}

// writerOpts turns a spec into tabix writer options, or reports an unknown preset.
func (t tabixSpec) writerOpts() (*tabix.WriterOpts, error) {
	opts := tabix.NewWriterOpts()
	switch strings.ToLower(t.preset) {
	case "bed":
		opts = opts.BED()
	case "vcf":
		opts = opts.VCF()
	case "gff", "gtf":
		opts = opts.GFF()
	case "":
		opts = opts.Columns(t.colSeq, t.colBeg, t.colEnd)
		if t.zeroBased {
			opts = opts.ZeroBased()
		}
	default:
		return nil, fmt.Errorf("unknown preset %q (use bed, vcf, or gff)", t.preset)
	}
	if t.meta != "" {
		opts = opts.Meta(t.meta[0])
	}
	if t.skip > 0 {
		opts = opts.Skip(t.skip)
	}
	return opts, nil
}
