package vcfcmd

import (
	"errors"
	"fmt"
	"regexp"
	"strings"

	"github.com/compgenlab/cghts/htsio/tabix"
	"github.com/compgenlab/cghts/vcf"
	"github.com/compgenlab/cgkit/internal/locator"
	"github.com/spf13/cobra"
)

// vcfTbi backs the --tbi flag shared by the VCF-writing commands. Only one
// subcommand runs per process, so a single variable is enough.
var vcfTbi bool

// addVcfOutputFlags registers the -o/--output and --tbi flags shared by the VCF
// commands that write a VCF stream.
func addVcfOutputFlags(cmd *cobra.Command, output *string) {
	f := cmd.Flags()
	f.StringVarP(output, "output", "o", "-", "Output filename (bgzip-compressed if it ends in .gz or .bgz; - for stdout)")
	f.BoolVar(&vcfTbi, "tbi", false, "Also write a tabix (.tbi) index for the output (requires -o with a .gz/.bgz name)")
}

// writeVcfTbi writes a companion .tbi index for a finished BGZF VCF file. Only
// vcf-concat and vcf-merge produce sorted output by construction; the rest pass
// their input's order through, so the index pass refuses an unsorted file rather
// than building an index whose linear offsets lie. The output itself is already
// written and stays valid; only the index is skipped.
func writeVcfTbi(filename string) error {
	err := tabix.NewIndexWriter(tabix.NewWriterOpts().VCF()).WriteIndex(filename)
	var unsorted *tabix.UnsortedError
	if errors.As(err, &unsorted) {
		return fmt.Errorf("wrote %s but skipped --tbi: %w\n"+
			"       sort and index it with: cgkit tab-sort -p vcf -o SORTED.vcf.gz %s", filename, err, filename)
	}
	return err
}

// openVcfWriter returns a VcfWriter for output, writing to stdout when output is
// "" or "-". The returned closer is nil for stdout (call writer.Close instead).
// With --tbi, the closer writes a companion .tbi index for the finished file.
func openVcfWriter(cmd *cobra.Command, output string) (*vcf.VcfWriter, func() error, error) {
	if output == "" || output == "-" {
		if vcfTbi {
			return nil, nil, fmt.Errorf("--tbi requires an output file (-o NAME.vcf.gz)")
		}
		return vcf.NewVcfWriter(cmd.OutOrStdout()), nil, nil
	}
	if vcfTbi && !strings.HasSuffix(output, ".gz") && !strings.HasSuffix(output, ".bgz") {
		return nil, nil, fmt.Errorf("--tbi requires a bgzip output name ending in .gz or .bgz, got %q", output)
	}
	if err := locator.CheckLocalOutput("-o/--output", output); err != nil {
		return nil, nil, err
	}
	w, err := vcf.OpenVcfWriter(output)
	if err != nil {
		return nil, nil, err
	}
	if !vcfTbi {
		return w, w.Close, nil
	}
	closeFn := func() error {
		if err := w.Close(); err != nil {
			return err
		}
		return writeVcfTbi(output)
	}
	return w, closeFn, nil
}

// globMatch reports whether s matches the glob pattern, porting ngsutilsj's
// GlobUtils: the pattern is anchored, '*' becomes ".*" and '?' becomes ".?".
func globMatch(s, glob string) bool {
	var b strings.Builder
	b.WriteByte('^')
	for _, c := range glob {
		switch c {
		case '*':
			b.WriteString(".*")
		case '?':
			b.WriteString(".?")
		case '.', '(', ')', '+', '|', '^', '$', '@', '%', '\\':
			b.WriteByte('\\')
			b.WriteRune(c)
		default:
			b.WriteRune(c)
		}
	}
	b.WriteByte('$')
	ok, err := regexp.MatchString(b.String(), s)
	return err == nil && ok
}
