package vcfcmd

import (
	"fmt"
	"io"
	"strings"

	"github.com/compgenlab/cghts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfRemoveFlagsOutput string
	vcfRemoveFlagsKey    string
	vcfRemoveFlagsAlways bool
)

var vcfRemoveFlagsCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-remove-flags <input.vcf>",
	Short:       "Replace all INFO flags with a comma-separated list",
	Long: `Replace every Flag-typed INFO field with a single key holding the set flags as
a comma-separated list (FOO;BAR => FLAGS=FOO,BAR).

  --key NAME   name for the new INFO key (default FLAGS)
  --always     always emit the key, using "." when no flags are set`,
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

		var flagIDs []string
		flagSet := map[string]bool{}
		for _, id := range header.InfoIDs() {
			if d, ok := header.InfoDef(id); ok && strings.EqualFold(d.Type, "Flag") {
				flagIDs = append(flagIDs, id)
				flagSet[id] = true
			}
		}
		if len(flagIDs) == 0 {
			return fmt.Errorf("no INFO flags defined in VCF")
		}
		for _, id := range flagIDs {
			header.RemoveInfo(id)
		}
		header.AddInfo(&vcf.AnnotationDef{
			IsInfo: true, ID: vcfRemoveFlagsKey, Number: ".", Type: "String",
			Description: "INFO Flag values as CSV (" + strings.Join(flagIDs, ",") + ")",
		})
		stampVcfProvenance(header, "vcf-remove-flags")

		writer, closeFn, err := openVcfWriter(cmd, vcfRemoveFlagsOutput)
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
			var setFlags []string
			for _, id := range rec.Info().Keys() {
				if flagSet[id] {
					setFlags = append(setFlags, id)
				}
			}
			if len(setFlags) > 0 {
				for _, id := range setFlags {
					rec.Info().Remove(id)
				}
				rec.AddInfo(vcfRemoveFlagsKey, strings.Join(setFlags, ","))
			} else if vcfRemoveFlagsAlways {
				rec.AddInfo(vcfRemoveFlagsKey, ".")
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
	f := vcfRemoveFlagsCmd.Flags()
	f.StringVarP(&vcfRemoveFlagsOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
	f.StringVar(&vcfRemoveFlagsKey, "key", "FLAGS", "Name for the new INFO key")
	f.BoolVar(&vcfRemoveFlagsAlways, "always", false, "Always include the key (use '.' when no flags set)")
}
