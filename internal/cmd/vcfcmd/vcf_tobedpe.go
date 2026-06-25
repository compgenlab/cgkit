package vcfcmd

import (
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/compgenlab/hts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfToBedpeOutput      string
	vcfToBedpePassing     bool
	vcfToBedpeNoDelOffset bool
	vcfToBedpeUniqueEvent bool
	vcfToBedpeAltChrom    string
	vcfToBedpeAltPos      string
	vcfToBedpeName        string
	vcfToBedpeScore       string
)

var vcfToBedpeCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-tobedpe <input.vcf>",
	Short:       "Convert a structural-variant VCF to BEDPE format",
	Long: `Convert a structural-variant VCF to BEDPE (paired-end BED). Each ALT breakpoint
becomes a line of two intervals: chrom1/start1/stop1 and chrom2/start2/stop2.

  --alt-chrom KEY   INFO field for the partner chromosome (default: from ALT)
  --alt-pos KEY     INFO field for the partner position (default: END / from ALT)
  --no-del-offset   don't offset deletion coordinates by one base
  --unique-event    one line per EVENT (requires an EVENT INFO field)
  --name KEY        INFO field (or @ID) to use as the BEDPE name
  --score KEY[:SAMPLE[:ALLELE]]   INFO/FORMAT field to use as the score
  --passing         only output passing variants`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		reader, err := openVcfInput(cmd, args[0])
		if err != nil {
			return err
		}
		defer reader.Close()
		header, err := reader.Header()
		if err != nil {
			return err
		}
		if vcfToBedpeUniqueEvent {
			if _, ok := header.InfoDef("EVENT"); !ok {
				return fmt.Errorf("--unique-event requires an EVENT INFO annotation")
			}
		}

		// Parse --score into KEY[:SAMPLE[:ALLELE]].
		var scoreKey, scoreSample, scoreAllele string
		scoreSampleIdx := -1
		if vcfToBedpeScore != "" {
			spl := strings.Split(vcfToBedpeScore, ":")
			scoreKey = spl[0]
			if len(spl) > 1 {
				scoreSample = spl[1]
				if scoreSample != "" && scoreSample != "INFO" {
					scoreSampleIdx = header.SampleIndex(scoreSample)
					if scoreSampleIdx < 0 {
						return fmt.Errorf("--score: sample not found: %s", scoreSample)
					}
				}
			}
			if len(spl) > 2 {
				scoreAllele = spl[2]
			}
		}

		out, closeFn, err := openOutput(cmd, vcfToBedpeOutput)
		if err != nil {
			return err
		}
		hdr := []string{"#chrom1", "start1", "stop1", "chrom2", "start2", "stop2", "name"}
		if vcfToBedpeScore != "" {
			hdr = append(hdr, "score")
		}
		fmt.Fprintln(out, strings.Join(hdr, "\t"))

		delOffset := !vcfToBedpeNoDelOffset
		seenEvents := map[string]bool{}
		for {
			rec, err := reader.NextRecord()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if vcfToBedpePassing && rec.IsFiltered() {
				continue
			}
			if vcfToBedpeUniqueEvent {
				if v, ok := rec.Info().Get("EVENT"); ok && v.String() != "" {
					if seenEvents[v.String()] {
						continue
					}
					seenEvents[v.String()] = true
				}
			}
			for _, alt := range rec.AltPositions(vcfToBedpeAltChrom, vcfToBedpeAltPos, "", "") {
				row := bedpeCoords(rec, alt, delOffset)
				row = append(row, bedpeName(rec, alt))
				if vcfToBedpeScore != "" {
					s, err := bedpeScore(rec, scoreKey, scoreSample, scoreSampleIdx, scoreAllele)
					if err != nil {
						return err
					}
					row = append(row, s)
				}
				fmt.Fprintln(out, strings.Join(row, "\t"))
			}
		}
		if closeFn != nil {
			return closeFn()
		}
		return nil
	},
}

// bedpeCoords builds the six coordinate columns for one breakpoint.
func bedpeCoords(rec *vcf.VcfRecord, alt vcf.AltPos, delOffset bool) []string {
	pos := rec.Pos
	var s1, e1, s2, e2 int
	switch {
	case delOffset && alt.Type == vcf.VarDEL:
		s1, e1 = pos, pos+1
		s2, e2 = alt.Pos-1, alt.Pos
	default: // INS and everything else use the same point/interval form
		s1, e1 = pos-1, pos
		s2, e2 = alt.Pos-1, alt.Pos
	}
	return []string{
		rec.Chrom, strconv.Itoa(s1), strconv.Itoa(e1),
		alt.Chrom, strconv.Itoa(s2), strconv.Itoa(e2),
	}
}

// bedpeName resolves the BEDPE name column.
func bedpeName(rec *vcf.VcfRecord, alt vcf.AltPos) string {
	switch {
	case vcfToBedpeName == "@ID":
		return idOrMissing(rec)
	case vcfToBedpeName != "":
		if v, ok := rec.Info().Get(vcfToBedpeName); ok {
			return v.String()
		}
		return "."
	case alt.Type == vcf.VarDEL:
		return "<DEL>"
	case alt.Type == vcf.VarINS:
		return "<INS>"
	default:
		return strings.Join(rec.Alt(), ",")
	}
}

// bedpeScore resolves the BEDPE score column from an INFO or FORMAT field.
func bedpeScore(rec *vcf.VcfRecord, key, sample string, sampleIdx int, allele string) (string, error) {
	if sample == "" || sample == "INFO" {
		if v, ok := rec.Info().Get(key); ok {
			return v.String(), nil
		}
		return ".", nil
	}
	s, err := rec.Sample(sampleIdx)
	if err != nil {
		return "", err
	}
	v, ok := s.Get(key)
	if !ok {
		return ".", nil
	}
	return v.StringFor(allele)
}

func init() {
	f := vcfToBedpeCmd.Flags()
	f.StringVarP(&vcfToBedpeOutput, "output", "o", "-", "Output filename (- for stdout)")
	f.BoolVar(&vcfToBedpePassing, "passing", false, "Only output passing variants")
	f.BoolVar(&vcfToBedpeNoDelOffset, "no-del-offset", false, "Don't offset deletion coordinates by one base")
	f.BoolVar(&vcfToBedpeUniqueEvent, "unique-event", false, "Only output one set of coordinates per EVENT")
	f.StringVar(&vcfToBedpeAltChrom, "alt-chrom", "", "Use an alternate INFO field for the partner chromosome")
	f.StringVar(&vcfToBedpeAltPos, "alt-pos", "END", "Use an alternate INFO field for the partner position")
	f.StringVar(&vcfToBedpeName, "name", "", "INFO field (or @ID) to use as the BEDPE name")
	f.StringVar(&vcfToBedpeScore, "score", "", "INFO/FORMAT field for the score (KEY[:SAMPLE[:ALLELE]])")
}
