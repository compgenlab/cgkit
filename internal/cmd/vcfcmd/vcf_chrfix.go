package vcfcmd

import (
	"fmt"
	"io"
	"strings"

	"github.com/compgenlab/hts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfChrFixOutput  string
	vcfChrFixUCSC    bool
	vcfChrFixEnsembl bool
	vcfChrFixPrimary bool
	vcfChrFixContigs string
)

var vcfChrFixCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-chrfix <input.vcf>",
	Short:       "Change the reference (chrom) format (Ensembl/UCSC)",
	Long: `Convert chromosome names between UCSC (chr1, chrM) and Ensembl (1, MT) forms,
and optionally drop non-primary contigs.

  --ucsc            convert to UCSC names (chr1, chrM)
  --ensembl         convert to Ensembl names (1, MT)
  --primary-human   keep only primary human contigs (1-22, X, Y, M)
  --contigs CSV     keep only these contigs (matched after conversion)`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if !vcfChrFixUCSC && !vcfChrFixEnsembl && !vcfChrFixPrimary && vcfChrFixContigs == "" {
			return fmt.Errorf("no changes specified (use --ucsc, --ensembl, --primary-human, or --contigs)")
		}
		if vcfChrFixUCSC && vcfChrFixEnsembl {
			return fmt.Errorf("you must set either --ucsc or --ensembl, not both")
		}

		var keepContigs map[string]bool
		if vcfChrFixContigs != "" {
			keepContigs = map[string]bool{}
			for _, c := range strings.Split(vcfChrFixContigs, ",") {
				keepContigs[strings.TrimSpace(c)] = true
			}
		}
		mapName := func(name string) string {
			switch {
			case vcfChrFixUCSC:
				return vcf.ToUCSC(name)
			case vcfChrFixEnsembl:
				return vcf.ToEnsembl(name)
			}
			return name
		}
		keepChrom := func(chrom string) bool {
			if keepContigs != nil {
				return keepContigs[chrom]
			}
			if vcfChrFixPrimary {
				return vcf.IsPrimaryHuman(chrom)
			}
			return true
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
		// Rename header contigs through the mapping.
		for _, old := range append([]string(nil), header.ContigNames()...) {
			newID := mapName(old)
			if newID == old {
				continue
			}
			length := int64(-1)
			if d, ok := header.ContigDef(old); ok {
				length = d.Length
			}
			header.RemoveContig(old)
			header.AddContig(&vcf.ContigDef{ID: newID, Length: length})
		}
		// Drop contigs that won't be kept.
		if vcfChrFixPrimary || keepContigs != nil {
			for _, c := range append([]string(nil), header.ContigNames()...) {
				if !keepChrom(c) {
					header.RemoveContig(c)
				}
			}
		}
		stampVcfProvenance(header, "vcf-chrfix")

		writer, closeFn, err := openVcfWriter(cmd, vcfChrFixOutput)
		if err != nil {
			return err
		}
		if err := writer.WriteHeader(header); err != nil {
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
			if newChrom := mapName(rec.Chrom); newChrom != rec.Chrom {
				rec.SetChrom(newChrom)
			}
			if !keepChrom(rec.Chrom) {
				continue
			}
			if err := writer.WriteRecord(rec); err != nil {
				return err
			}
		}
		if closeFn != nil {
			return closeFn()
		}
		return writer.Close()
	},
}

func init() {
	f := vcfChrFixCmd.Flags()
	f.StringVarP(&vcfChrFixOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
	f.BoolVar(&vcfChrFixUCSC, "ucsc", false, "Convert to UCSC references (chr1, chr2, ...)")
	f.BoolVar(&vcfChrFixEnsembl, "ensembl", false, "Convert to Ensembl references (1, 2, ...)")
	f.BoolVar(&vcfChrFixPrimary, "primary-human", false, "Keep only primary human contigs (1-22, X, Y, M)")
	f.StringVar(&vcfChrFixContigs, "contigs", "", "Keep only these contigs (comma-separated)")
}
