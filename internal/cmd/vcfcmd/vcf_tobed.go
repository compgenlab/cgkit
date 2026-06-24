package vcfcmd

import (
	"bufio"
	"fmt"
	"io"

	"github.com/spf13/cobra"
)

var (
	vcfToBedOutput     string
	vcfToBedPassing    bool
	vcfToBedIncludePos bool
	vcfToBedPadding    int
	vcfToBedAltChrom   string
	vcfToBedAltPos     string
	vcfToBedRegion     string
)

var vcfToBedCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-tobed <input.vcf>",
	Short:       "Export allele positions from a VCF file to BED format",
	Long: `Export allele positions from a VCF file to BED format.

Each variant is written as a BED interval [POS-1, end), where end is the
alternate-allele position (for SVs this is resolved from the ALT field or an
INFO key such as END). The fourth column is the variant type (SNV, DEL, BND,
...), or the original CHROM_POS when --include-pos is given.

Breakends (BND) that span two different chromosomes cannot be represented in
BED and are skipped.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		src, err := openRecordSource(cmd, args[0], vcfToBedRegion)
		if err != nil {
			return err
		}
		defer src.close()

		w, closeFn, err := openOutput(cmd, vcfToBedOutput)
		if err != nil {
			return err
		}
		out := bufio.NewWriter(w)

		for {
			rec, err := src.next()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if vcfToBedPassing && rec.IsFiltered() {
				continue
			}
			chrom := rec.Chrom
			pos := rec.Pos
			for _, alt := range rec.AltPositions(vcfToBedAltChrom, vcfToBedAltPos, "", "") {
				chrom2 := alt.Chrom
				if vcfToBedAltChrom != "" {
					if v, ok := rec.Info().Get(vcfToBedAltChrom); ok {
						chrom2 = v.String()
					}
				}
				if chrom2 != chrom {
					// BND across chromosomes cannot be written to BED.
					continue
				}
				endpos := alt.Pos
				if vcfToBedAltPos != "" {
					if v, ok := rec.Info().Get(vcfToBedAltPos); ok {
						if n, err := v.Int(); err == nil {
							endpos = n
						}
					}
				}
				name := alt.Type.String()
				if vcfToBedIncludePos {
					name = fmt.Sprintf("%s_%d", rec.Chrom, rec.Pos)
				}
				fmt.Fprintf(out, "%s\t%d\t%d\t%s\n", chrom, (pos-1)-vcfToBedPadding, endpos+vcfToBedPadding, name)
			}
		}

		if err := out.Flush(); err != nil {
			return err
		}
		if closeFn != nil {
			return closeFn()
		}
		return nil
	},
}

func init() {
	vcfToBedCmd.Flags().StringVarP(&vcfToBedOutput, "output", "o", "-", "Output filename (- for stdout)")
	vcfToBedCmd.Flags().BoolVar(&vcfToBedPassing, "passing", false, "Only output passing variants")
	vcfToBedCmd.Flags().BoolVar(&vcfToBedIncludePos, "include-pos", false, "Use CHROM_POS as the name field (without padding)")
	vcfToBedCmd.Flags().IntVar(&vcfToBedPadding, "padding", 0, "Add extra padding on either side")
	vcfToBedCmd.Flags().StringVar(&vcfToBedAltChrom, "alt-chrom", "", "Use an alternate INFO field for the chromosome (default: extracted from ALT)")
	vcfToBedCmd.Flags().StringVar(&vcfToBedAltPos, "alt-pos", "", "Use an alternate INFO field for the position (default: extracted from ALT, or END)")
	vcfToBedCmd.Flags().StringVar(&vcfToBedRegion, "region", "", "Only variants in this 1-based region (chrom:start-end, or chrom); requires a tabix-indexed file")
}
