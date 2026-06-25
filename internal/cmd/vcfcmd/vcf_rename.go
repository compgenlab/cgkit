package vcfcmd

import (
	"fmt"
	"io"
	"strings"

	"github.com/spf13/cobra"
)

var (
	vcfRenameOutput  string
	vcfRenameSamples []string
	vcfRenameFiles   []string
)

var vcfRenameCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-rename <input.vcf>",
	Short:       "Change the names of samples",
	Long: `Rename samples in a VCF file. The sample columns and all record data are
unchanged; only the sample names in the header are updated.

  --sample OLD:NEW    rename a sample (OLD may be a 1-based number; repeatable)
  --samples-file, -f  tab-delimited file of "oldname<TAB>newname" lines`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		var oldNames, newNames []string
		for _, file := range vcfRenameFiles {
			lines, err := readLines(file)
			if err != nil {
				return err
			}
			for _, line := range lines {
				cols := strings.Split(line, "\t")
				if len(cols) == 2 {
					oldNames = append(oldNames, cols[0])
					newNames = append(newNames, cols[1])
				}
			}
		}
		for _, v := range vcfRenameSamples {
			spl := strings.Split(v, ":")
			if len(spl) != 2 {
				return fmt.Errorf("invalid --sample %q (want OLD:NEW)", v)
			}
			oldNames = append(oldNames, spl[0])
			newNames = append(newNames, spl[1])
		}
		if len(oldNames) == 0 {
			return fmt.Errorf("you must specify at least one sample to rename")
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
		for i := range oldNames {
			if err := header.RenameSample(oldNames[i], newNames[i]); err != nil {
				return err
			}
		}
		stampVcfProvenance(header, "vcf-rename")

		writer, closeFn, err := openVcfWriter(cmd, vcfRenameOutput)
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
	f := vcfRenameCmd.Flags()
	f.StringVarP(&vcfRenameOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
	f.StringArrayVar(&vcfRenameSamples, "sample", nil, "Rename a sample as OLD:NEW (OLD may be a 1-based number; repeatable)")
	f.StringArrayVarP(&vcfRenameFiles, "samples-file", "f", nil, "Tab-delimited file of oldname<TAB>newname (repeatable)")
}
