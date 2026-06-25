package vcfcmd

import (
	"fmt"
	"strconv"

	"github.com/spf13/cobra"
)

var (
	vcfHeaderInfoInfo    bool
	vcfHeaderInfoFormat  bool
	vcfHeaderInfoSample  bool
	vcfHeaderInfoFilters bool
	vcfHeaderInfoContig  bool
)

var vcfHeaderInfoCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-header-info <input.vcf>",
	Short:       "Extract annotation/named fields from a VCF header",
	Long: `Write header metadata as tab-delimited text. Specify exactly one field type:

  --info / --format / --filters   id<TAB>description (one per line)
  --sample                        sample name (one per line)
  --contig                        id<TAB>length (one per line)`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		n := 0
		for _, b := range []bool{vcfHeaderInfoInfo, vcfHeaderInfoFormat, vcfHeaderInfoSample, vcfHeaderInfoFilters, vcfHeaderInfoContig} {
			if b {
				n++
			}
		}
		if n != 1 {
			return fmt.Errorf("specify exactly one of --info, --format, --sample, --filters, --contig")
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
		out := cmd.OutOrStdout()
		switch {
		case vcfHeaderInfoInfo:
			for _, id := range header.InfoIDs() {
				d, _ := header.InfoDef(id)
				fmt.Fprintf(out, "%s\t%s\n", d.ID, d.Description)
			}
		case vcfHeaderInfoFormat:
			for _, id := range header.FormatIDs() {
				d, _ := header.FormatDef(id)
				fmt.Fprintf(out, "%s\t%s\n", d.ID, d.Description)
			}
		case vcfHeaderInfoSample:
			for _, s := range header.Samples() {
				fmt.Fprintln(out, s)
			}
		case vcfHeaderInfoFilters:
			for _, id := range header.FilterIDs() {
				d, _ := header.FilterDef(id)
				fmt.Fprintf(out, "%s\t%s\n", d.ID, d.Description)
			}
		case vcfHeaderInfoContig:
			for _, id := range header.ContigNames() {
				d, _ := header.ContigDef(id)
				length := ""
				if d.Length >= 0 {
					length = strconv.FormatInt(d.Length, 10)
				}
				fmt.Fprintf(out, "%s\t%s\n", d.ID, length)
			}
		}
		return nil
	},
}

func init() {
	f := vcfHeaderInfoCmd.Flags()
	f.BoolVar(&vcfHeaderInfoInfo, "info", false, "Output INFO field definitions")
	f.BoolVar(&vcfHeaderInfoFormat, "format", false, "Output FORMAT field definitions")
	f.BoolVar(&vcfHeaderInfoSample, "sample", false, "Output sample names")
	f.BoolVar(&vcfHeaderInfoFilters, "filters", false, "Output FILTER definitions")
	f.BoolVar(&vcfHeaderInfoContig, "contig", false, "Output contig definitions")
}
