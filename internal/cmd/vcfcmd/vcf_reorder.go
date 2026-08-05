package vcfcmd

import (
	"bufio"
	"fmt"
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
1-based number. Samples omitted from the new order are dropped.

A requested sample the file does not have is an error, and the message lists
the samples it does have. This used to warn and carry on, which quietly
produced a VCF one column short of the cohort that was asked for.

FORMAT values are not parsed: the sample columns are moved verbatim.`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
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

		var order []int
		return runVcfStream(cmd, vcfStream{
			name: "vcf-reorder",
			in:   args[0],
			out:  vcfReorderOutput,
			header: func(header *vcf.VcfHeader) error {
				orig := header.Samples()

				// An unresolvable name is fatal. It used to warn and carry on, which is
				// the worst of the options: naming one sample wrong silently produced a
				// VCF one column short of the cohort that was asked for, and naming them
				// all wrong produced a header with no FORMAT and no sample columns over
				// records that still had theirs -- a file that is not valid VCF at all,
				// written with exit status 0.
				var newNames []string
				var missing []string
				for _, name := range requested {
					idx := header.SampleIndex(name)
					if idx < 0 || idx >= len(orig) {
						missing = append(missing, name)
						continue
					}
					order = append(order, idx)
					newNames = append(newNames, name)
				}
				if len(missing) > 0 {
					return fmt.Errorf("no such sample%s: %s\n  this file has: %s",
						plural(len(missing)), strings.Join(missing, ", "), strings.Join(orig, ", "))
				}

				header.SetSamples(newNames)
				return nil
			},
			record: func(*vcf.VcfRecord) (bool, error) { return true, nil },
			write: func(w *vcf.VcfWriter, rec *vcf.VcfRecord) error {
				return w.WriteLine(rec.ReorderSamplesLine(order))
			},
		})
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
	addVcfOutputFlags(vcfReorderCmd, &vcfReorderOutput)
	vcfReorderCmd.Flags().StringArrayVarP(&vcfReorderSamples, "sample", "s", nil, "New sample order (comma-separated, repeatable)")
	vcfReorderCmd.Flags().StringVar(&vcfReorderSamplesFile, "samples-file", "", "File with the new sample order, one per line")
}
