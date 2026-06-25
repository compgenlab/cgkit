package vcfcmd

import (
	"fmt"
	"io"
	"strconv"

	"github.com/compgenlab/hts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfSplitOut string
	vcfSplitNum int
)

var vcfSplitCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-split <input.vcf>",
	Short:       "Split a VCF file into smaller files with N variants each",
	Long: `Split a VCF file into multiple bgzipped files of N variants each. Each output
file gets a fresh copy of the header. Outputs are named BASE.1.vcf.gz,
BASE.2.vcf.gz, and so on.

  --out BASE    base output name (required)
  --num N       variants per output file (required)`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if vcfSplitOut == "" {
			return fmt.Errorf("you must specify a base output name with --out")
		}
		if vcfSplitNum <= 0 {
			return fmt.Errorf("--num must be a positive number")
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
		stampVcfProvenance(header, "vcf-split")

		var writer *vcf.VcfWriter
		fileNum, inFile := 0, 0
		closeChunk := func() error {
			if writer == nil {
				return nil
			}
			err := writer.Close()
			writer = nil
			return err
		}

		for {
			rec, err := reader.NextRecord()
			if err == io.EOF {
				break
			}
			if err != nil {
				closeChunk()
				return err
			}
			if writer == nil {
				fileNum++
				w, oerr := vcf.OpenVcfWriter(vcfSplitOut + "." + strconv.Itoa(fileNum) + ".vcf.gz")
				if oerr != nil {
					return oerr
				}
				writer = w
				if herr := writer.WriteHeader(header); herr != nil {
					closeChunk()
					return herr
				}
				inFile = 0
			}
			if err := writer.WriteRecord(rec); err != nil {
				closeChunk()
				return err
			}
			inFile++
			if inFile >= vcfSplitNum {
				if err := closeChunk(); err != nil {
					return err
				}
			}
		}
		return closeChunk()
	},
}

func init() {
	f := vcfSplitCmd.Flags()
	f.StringVar(&vcfSplitOut, "out", "", "Base output name (outputs are BASE.N.vcf.gz)")
	f.IntVar(&vcfSplitNum, "num", 0, "Number of variants per output file")
}
