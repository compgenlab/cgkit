package vcfcmd

import (
	"bufio"
	"fmt"
	"io"

	"github.com/spf13/cobra"
)

var (
	vcfTsTvPassing bool
	vcfTsTvRegion  string
)

var vcfTsTvCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-tstv <input.vcf>",
	Short:       "Calculate a Ts/Tv ratio for SNVs",
	Long:        "Calculate the transition/transversion (Ts/Tv) ratio for SNVs in a VCF file.",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		src, err := openRecordSource(cmd, args[0], vcfTsTvRegion)
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
			if vcfTsTvPassing && rec.IsFiltered() {
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
// does, including "Infinity" and "NaN" for division by zero.
func javaRatio(num, den float64) string {
	if den == 0 {
		if num == 0 {
			return "NaN"
		}
		return "Infinity"
	}
	return javaDouble(num / den)
}

func init() {
	vcfTsTvCmd.Flags().BoolVar(&vcfTsTvPassing, "passing", false, "Only use passing variants")
	vcfTsTvCmd.Flags().StringVar(&vcfTsTvRegion, "region", "", "Only variants in this 1-based region (chrom:start-end, or chrom); requires a tabix-indexed file")
}
