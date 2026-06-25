package vcfcmd

import (
	"regexp"
	"strings"

	"github.com/compgenlab/hts/vcf"
	"github.com/spf13/cobra"
)

// openVcfWriter returns a VcfWriter for output, writing to stdout when output is
// "" or "-". The returned closer is nil for stdout (call writer.Close instead).
func openVcfWriter(cmd *cobra.Command, output string) (*vcf.VcfWriter, func() error, error) {
	if output == "" || output == "-" {
		return vcf.NewVcfWriter(cmd.OutOrStdout()), nil, nil
	}
	w, err := vcf.OpenVcfWriter(output)
	if err != nil {
		return nil, nil, err
	}
	return w, w.Close, nil
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
