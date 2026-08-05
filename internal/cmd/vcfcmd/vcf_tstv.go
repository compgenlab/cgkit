package vcfcmd

import (
	"bufio"
	"fmt"
	"io"

	"github.com/spf13/cobra"
)

var vcfTsTvCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-tstv <input.vcf>",
	Short:       "Calculate a Ts/Tv ratio for SNVs",
	Long:        "Calculate the transition/transversion (Ts/Tv) ratio for SNVs in a VCF file.",
	Args:        cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		src, err := openRecordSource(cmd, args[0], vcfRegion)
		if err != nil {
			return err
		}
		defer src.close()

		var tsCount, tvCount int64
		for {
			rec, err := src.next()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if vcfPassing && rec.IsFiltered() {
				continue
			}
			switch rec.CalcTsTv() {
			case -1:
				tsCount++
			case 1:
				tvCount++
			}
		}

		out := bufio.NewWriter(cmd.OutOrStdout())
		fmt.Fprintf(out, "Transitions (Ts)\t%d\n", tsCount)
		fmt.Fprintf(out, "Transversions (Tv)\t%d\n", tvCount)
		fmt.Fprintf(out, "Ts/Tv ratio\t%s\n", javaRatio(float64(tsCount), float64(tvCount)))
		return out.Flush()
	},
}

// javaRatio formats num/den the way Java's string concatenation of a double
// does, which is what ngsutilsj emits and what the parity tests compare against.
//
// The zero-denominator case deliberately diverges. Java gives "Infinity" or
// "NaN", and vcf-stats gave an empty string for the same situation -- one of
// which reads as a real number and the other as a bug, when what actually
// happened is that there were no transversions to divide by. Both report "-"
// now, which is the same thing the tabular commands use for a value that does
// not exist. testdata/sample.vcf has no transversions at all, so the parity
// tests do cover this -- they normalize the reference output for this one field
// rather than pretending the divergence is not there.
func javaRatio(num, den float64) string {
	if den == 0 {
		return "-"
	}
	return javaDouble(num / den)
}

func init() {
	addPassingFlag(vcfTsTvCmd, "Only use passing variants")
	addRegionFlag(vcfTsTvCmd)
}
