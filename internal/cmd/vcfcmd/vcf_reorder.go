package vcfcmd

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/compgenlab/cghts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfReorderOutput      string
	vcfReorderSamples     []string
	vcfReorderSamplesFile string
)

var vcfReorderCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-reorder <input.vcf>",
	Short:       "Reorder (or subset) the samples in a VCF file",
	Long: `Reorder the samples in a VCF file.

The new sample order is given with --sample (repeatable, comma-separated) or
--samples-file (one sample per line). Samples may be named or referenced by
1-based number. Samples omitted from the new order are dropped; a requested
sample that is not present is skipped with a warning.

FORMAT values are not parsed: the sample columns are moved verbatim.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		haveList := len(vcfReorderSamples) > 0
		haveFile := vcfReorderSamplesFile != ""
		if haveList == haveFile {
			return fmt.Errorf("you must specify exactly one of --sample or --samples-file")
		}

		var requested []string
		if haveFile {
			lines, err := readLines(vcfReorderSamplesFile)
			if err != nil {
				return err
			}
			requested = lines
		} else {
			for _, val := range vcfReorderSamples {
				for _, s := range strings.Split(val, ",") {
					requested = append(requested, strings.TrimSpace(s))
				}
			}
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
		orig := header.Samples()

		var order []int
		var newNames []string
		for _, name := range requested {
			idx := header.SampleIndex(name)
			if idx < 0 || idx >= len(orig) {
				fmt.Fprintf(cmd.ErrOrStderr(), "Missing sample: %s\n", name)
				continue
			}
			order = append(order, idx)
			newNames = append(newNames, name)
		}

		header.SetSamples(newNames)
		stampVcfProvenance(header, "vcf-reorder")

		var writer *vcf.VcfWriter
		var closeErr func() error
		if vcfReorderOutput == "" || vcfReorderOutput == "-" {
			writer = vcf.NewVcfWriter(cmd.OutOrStdout())
		} else {
			w, err := vcf.OpenVcfWriter(vcfReorderOutput)
			if err != nil {
				return err
			}
			writer = w
			closeErr = w.Close
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
			if err := writer.WriteLine(rec.ReorderSamplesLine(order)); err != nil {
				return err
			}
		}
		if closeErr != nil {
			return closeErr()
		}
		return writer.Close()
	},
}

func readLines(filename string) ([]string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var lines []string
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line != "" {
			lines = append(lines, line)
		}
	}
	return lines, sc.Err()
}

func init() {
	vcfReorderCmd.Flags().StringVarP(&vcfReorderOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
	vcfReorderCmd.Flags().StringArrayVarP(&vcfReorderSamples, "sample", "s", nil, "New sample order (comma-separated, repeatable)")
	vcfReorderCmd.Flags().StringVar(&vcfReorderSamplesFile, "samples-file", "", "File with the new sample order, one per line")
}
