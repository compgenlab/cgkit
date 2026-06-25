package vcfcmd

import (
	"fmt"
	"io"
	"strings"

	"github.com/compgenlab/hts/vcf"
	"github.com/spf13/cobra"
)

var vcfMergeOutput string

var vcfMergeCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-merge <input1.vcf> <input2.vcf> [input3.vcf ...]",
	Short:       "Combine VCF files with the same variants but different annotations",
	Long: `Merge VCF files that contain the same variants (same CHROM/POS/REF/ALT, in the
same order) but different annotations — for example a base VCF annotated by
several tools in parallel. The annotations are combined per variant; on a
conflict the first file on the command line wins. A variant missing from any
input is an error.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if len(args) < 2 {
			return fmt.Errorf("vcf-merge needs at least two input VCF files")
		}

		var streams []*recordSource
		closeAll := func() {
			for _, s := range streams {
				if s != nil {
					s.close()
				}
			}
		}
		for _, fn := range args {
			s, err := openConcatStream(cmd, fn, false)
			if err != nil {
				closeAll()
				return err
			}
			streams = append(streams, s)
		}
		defer closeAll()

		header := streams[0].header
		for _, s := range streams[1:] {
			if err := unionHeaderInto(header, s.header, true); err != nil {
				return err
			}
		}
		stampVcfProvenance(header, "vcf-merge")

		writer, closeFn, err := openVcfWriter(cmd, vcfMergeOutput)
		if err != nil {
			return err
		}
		if err := writer.WriteHeader(header); err != nil {
			return err
		}

		primary, secondaries := streams[0], streams[1:]
		for {
			rec, err := primary.next()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			for si, sec := range secondaries {
				next, err := sec.next()
				if err == io.EOF {
					return fmt.Errorf("%s is missing variant %s:%d", args[si+1], rec.Chrom, rec.Pos)
				}
				if err != nil {
					return err
				}
				if err := mergeRecord(rec, next, args[si+1]); err != nil {
					return err
				}
			}
			if err := writer.WriteRecord(rec); err != nil {
				return err
			}
		}
		// Every secondary must be exhausted too.
		for si, sec := range secondaries {
			if _, err := sec.next(); err != io.EOF {
				return fmt.Errorf("%s has more variants than the primary input", args[si+1])
			}
		}

		if closeFn != nil {
			return closeFn()
		}
		return writer.Close()
	},
}

// mergeRecord folds next's annotations into rec (primary wins on conflict),
// after verifying they describe the same variant.
func mergeRecord(rec, next *vcf.VcfRecord, srcName string) error {
	altR, altN := rec.Alt(), next.Alt()
	mismatch := rec.Chrom != next.Chrom || rec.Pos != next.Pos || rec.Ref != next.Ref || len(altR) != len(altN)
	if !mismatch {
		for i := range altR {
			if altR[i] != altN[i] {
				mismatch = true
				break
			}
		}
	}
	if mismatch {
		return fmt.Errorf("variants out of order: expected %s:%d %s>%s, %s has %s:%d %s>%s",
			rec.Chrom, rec.Pos, rec.Ref, strings.Join(altR, ","),
			srcName, next.Chrom, next.Pos, next.Ref, strings.Join(altN, ","))
	}

	// ID: first non-empty wins.
	if rec.ID() == "" && next.ID() != "" {
		rec.SetID(next.ID())
	}

	// FILTER: distinct union, primary codes first.
	present := map[string]bool{}
	for _, f := range rec.Filters() {
		present[f] = true
	}
	for _, f := range next.Filters() {
		if !present[f] {
			rec.AddFilter(f)
			present[f] = true
		}
	}

	// INFO: add keys absent from the primary.
	for _, k := range next.Info().Keys() {
		if !rec.Info().Contains(k) {
			v, _ := next.Info().Get(k)
			rec.AddInfo(k, v.String())
		}
	}

	// FORMAT: per sample, add keys absent from the primary.
	for i := 0; i < rec.NumSamples(); i++ {
		rs, err := rec.Sample(i)
		if err != nil {
			return err
		}
		ns, err := next.Sample(i)
		if err != nil {
			return err
		}
		for _, k := range ns.Keys() {
			if !rs.Contains(k) {
				v, _ := ns.Get(k)
				if err := rec.AddFormat(i, k, v.String()); err != nil {
					return err
				}
			}
		}
	}
	return nil
}

func init() {
	vcfMergeCmd.Flags().StringVarP(&vcfMergeOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
}
