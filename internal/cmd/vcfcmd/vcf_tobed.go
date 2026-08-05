package vcfcmd

import (
	"bufio"
	"fmt"
	"io"

	"github.com/spf13/cobra"

	"github.com/compgenlab/cghts/vcf"
)

var (
	vcfToBedOutput     string
	vcfToBedIncludePos bool
	vcfToBedPadding    int
	vcfToBedAltChrom   string
	vcfToBedAltPos     string
)

var vcfToBedCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-tobed <input.vcf>",
	Short:       "Export allele positions from a VCF file to BED format",
	Long: `Export allele positions from a VCF file to BED format.

Each variant is written as a BED interval [POS-1, end) covering the reference
bases the record spans: len(REF) for a plain variant, widened by INFO/END,
INFO/SVLEN or FORMAT/LEN where those apply. The fourth column is the variant
type (SNV, DEL, BND, ...), or the original CHROM_POS when --include-pos is
given.

gVCF reference blocks are skipped -- they describe coverage, not variants.

Breakends (BND) that span two different chromosomes cannot be represented in
BED and are skipped.`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		src, err := openRecordSource(cmd, args[0], vcfRegion)
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
			if vcfPassing && rec.IsFiltered() {
				continue
			}
			// A gVCF reference block is not a variant, so it has no place in a BED
			// of variant positions. Worth skipping explicitly: such a record used
			// to be emitted mislabelled as an SNV, or silently dropped when its ALT
			// was a bare ".".
			if rec.IsRefBlock() {
				continue
			}
			chrom := rec.Chrom
			pos := rec.Pos
			// The reference bases this record actually covers. AltPositions cannot
			// answer this: it resolves an SV's *partner breakpoint*, which may be on
			// another chromosome, where this asks how far the record itself reaches.
			// Using it for the end gave a plain deletion one base too many and
			// collapsed a multi-base-ALT deletion to a single base.
			_, refEnd := rec.RefSpan()
			for _, alt := range rec.AltPositions(vcfToBedAltChrom, vcfToBedAltPos, "", "") {
				// A gVCF variant record carries the block allele alongside the real
				// one ("G,<NON_REF>"), so the record itself is kept while this
				// allele is not -- it names no variant to report.
				if vcf.IsRefBlockAlt(alt.Alt) {
					continue
				}
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
				start, endpos := pos-1, refEnd
				// An explicit --alt-pos still wins: the caller is naming the INFO
				// field that holds the end, and for a breakend that partner may lie
				// upstream, so order the two rather than emit a negative interval.
				if vcfToBedAltPos != "" {
					if v, ok := rec.Info().Get(vcfToBedAltPos); ok {
						if n, err := v.Int(); err == nil {
							endpos = n
						}
					}
				} else if alt.Type == vcf.VarBND {
					endpos = alt.Pos
				}
				if endpos < start {
					start, endpos = endpos, start
				}
				name := alt.Type.String()
				if vcfToBedIncludePos {
					name = fmt.Sprintf("%s_%d", rec.Chrom, rec.Pos)
				}
				fmt.Fprintf(out, "%s\t%d\t%d\t%s\n", chrom, start-vcfToBedPadding, endpos+vcfToBedPadding, name)
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
	addOutputFlag(vcfToBedCmd, &vcfToBedOutput)
	addPassingFlag(vcfToBedCmd, "Only output passing variants")
	vcfToBedCmd.Flags().BoolVar(&vcfToBedIncludePos, "include-pos", false, "Use CHROM_POS as the name field (without padding)")
	vcfToBedCmd.Flags().IntVar(&vcfToBedPadding, "padding", 0, "Add extra padding on either side")
	vcfToBedCmd.Flags().StringVar(&vcfToBedAltChrom, "alt-chrom", "", "Use an alternate INFO field for the chromosome (default: extracted from ALT)")
	vcfToBedCmd.Flags().StringVar(&vcfToBedAltPos, "alt-pos", "", "Use an alternate INFO field for the position (default: extracted from ALT, or END)")
	addRegionFlag(vcfToBedCmd)
}
