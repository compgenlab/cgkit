package vcfcmd

import (
	"fmt"
	"io"
	"os"

	"github.com/compgenlab/hts/seqio"
	"github.com/compgenlab/hts/support/sequtils"
	"github.com/compgenlab/hts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfSvToFastaOutput     string
	vcfSvToFastaPassing    bool
	vcfSvToFastaIncludeRef bool
	vcfSvToFastaBND        bool
	vcfSvToFastaFlanking   int
	vcfSvToFastaSVType     string
	vcfSvToFastaCT         string
	vcfSvToFastaAltChrom   string
	vcfSvToFastaAltPos     string
)

var vcfSvToFastaCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-svtofasta <genome.fa> <input.vcf>",
	Short:       "Extract SV breakend flanking sequences to FASTA",
	Long: `For breakend (BND) structural variants, extract the flanking reference sequence
on each side of the breakpoint, join them in breakend orientation, and write the
result as FASTA. Requires an indexed reference (genome.fa.fai).

  --bnd             export breakend/translocation sequences (required for output)
  --flanking N      flanking bases to include on each side (default 1000)
  --include-ref     also write the wild-type reference sequences
  --svtype KEY      INFO field for the SV type (default SVTYPE)
  --ct KEY          INFO field for the connection type (5to5/5to3/3to3/3to5)
  --alt-chrom KEY   INFO field for the partner chromosome (default: from ALT)
  --alt-pos KEY     INFO field for the partner position (default END)
  --passing         only process passing variants`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if len(args) < 2 {
			return fmt.Errorf("vcf-svtofasta needs a genome FASTA and an input VCF")
		}
		fastaFile, vcfFile := args[0], args[1]
		if _, err := os.Stat(fastaFile + ".fai"); err != nil {
			return fmt.Errorf("missing FAI index for %s (run samtools faidx)", fastaFile)
		}

		ref, err := seqio.OpenReference(fastaFile)
		if err != nil {
			return err
		}
		defer ref.Close()

		reader, err := openVcfInput(cmd, vcfFile)
		if err != nil {
			return err
		}
		defer reader.Close()
		if _, err := reader.Header(); err != nil {
			return err
		}

		out, closeFn, err := openOutput(cmd, vcfSvToFastaOutput)
		if err != nil {
			return err
		}
		for {
			rec, err := reader.NextRecord()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if vcfSvToFastaPassing && rec.IsFiltered() {
				continue
			}
			for _, alt := range rec.AltPositions(vcfSvToFastaAltChrom, vcfSvToFastaAltPos, vcfSvToFastaSVType, vcfSvToFastaCT) {
				if vcfSvToFastaBND && alt.Type == vcf.VarBND {
					if err := writeBreakendFasta(out, ref, rec, alt); err != nil {
						return err
					}
				}
			}
		}
		if closeFn != nil {
			return closeFn()
		}
		return nil
	},
}

// writeBreakendFasta writes the joined flanking sequence for one BND breakpoint
// (and, with --include-ref, the wild-type references). It ports
// VCFSVToFASTA.writeBND.
func writeBreakendFasta(out io.Writer, ref seqio.ReferenceReader, rec *vcf.VcfRecord, alt vcf.AltPos) error {
	chrom, pos := rec.Chrom, rec.Pos
	chrom2, pos2 := alt.Chrom, alt.Pos
	flank := vcfSvToFastaFlanking

	fetch := func(name string, start, end int) (string, error) {
		b, err := ref.GetSequenceRange(name, start, end)
		return string(b), err
	}

	var seq string
	switch alt.ConnType {
	case vcf.Conn3to5:
		a, err := fetch(chrom, pos-flank, pos)
		if err != nil {
			return err
		}
		b, err := fetch(chrom2, pos2, pos2+flank)
		if err != nil {
			return err
		}
		seq = a + b
	case vcf.Conn5to3:
		a, err := fetch(chrom, pos, pos+flank)
		if err != nil {
			return err
		}
		b, err := fetch(chrom2, pos2-flank, pos2)
		if err != nil {
			return err
		}
		seq = b + a
	case vcf.Conn5to5:
		a, err := fetch(chrom, pos, pos+flank)
		if err != nil {
			return err
		}
		b, err := fetch(chrom2, pos2, pos2+flank)
		if err != nil {
			return err
		}
		seq = sequtils.ReverseComplement(b) + a
	case vcf.Conn3to3:
		a, err := fetch(chrom, pos-flank, pos)
		if err != nil {
			return err
		}
		b, err := fetch(chrom2, pos2-flank, pos2)
		if err != nil {
			return err
		}
		seq = a + sequtils.ReverseComplement(b)
	default:
		return nil
	}

	id := rec.ID()
	fmt.Fprintf(out, ">sv|%s|%d|%s|%d|%s|%s|%s\n%s\n",
		chrom, pos, chrom2, pos2, id, alt.Type, alt.ConnType, seq)

	if !vcfSvToFastaIncludeRef {
		return nil
	}
	refA, err := fetch(chrom, pos-flank, pos+flank)
	if err != nil {
		return err
	}
	fmt.Fprintf(out, ">ref|%s|%d|%d|%s|%s|%s|A\n%s\n",
		chrom, pos-flank, pos+flank, id, alt.Type, alt.ConnType, refA)

	if alt.ConnType == vcf.Conn5to3 || alt.ConnType == vcf.Conn3to5 {
		// Faithful to ngsutilsj, which fetches from chrom (not chrom2) here.
		refB, err := fetch(chrom, pos2-flank, pos2+flank)
		if err != nil {
			return err
		}
		fmt.Fprintf(out, ">ref|%s|%d|%d|%s|%s|%s|B\n%s\n",
			chrom2, pos2-flank, pos2+flank, id, alt.Type, alt.ConnType, refB)
	} else {
		refB, err := fetch(chrom2, pos2-flank, pos2+flank)
		if err != nil {
			return err
		}
		fmt.Fprintf(out, ">ref|%s|%d|%d|%s|%s|%s|B\n%s\n",
			chrom2, pos2+flank, pos2-flank, id, alt.Type, alt.ConnType, sequtils.ReverseComplement(refB))
	}
	return nil
}

func init() {
	f := vcfSvToFastaCmd.Flags()
	f.StringVarP(&vcfSvToFastaOutput, "output", "o", "-", "Output filename (- for stdout)")
	f.BoolVar(&vcfSvToFastaBND, "bnd", false, "Export breakend/translocation sequences (required for output)")
	f.IntVar(&vcfSvToFastaFlanking, "flanking", 1000, "Flanking bases to include on each side")
	f.BoolVar(&vcfSvToFastaIncludeRef, "include-ref", false, "Also write wild-type reference sequences")
	f.BoolVar(&vcfSvToFastaPassing, "passing", false, "Only process passing variants")
	f.StringVar(&vcfSvToFastaSVType, "svtype", "SVTYPE", "INFO field for the SV type")
	f.StringVar(&vcfSvToFastaCT, "ct", "", "INFO field for the connection type")
	f.StringVar(&vcfSvToFastaAltChrom, "alt-chrom", "", "Use an alternate INFO field for the partner chromosome")
	f.StringVar(&vcfSvToFastaAltPos, "alt-pos", "END", "Use an alternate INFO field for the partner position")
}
