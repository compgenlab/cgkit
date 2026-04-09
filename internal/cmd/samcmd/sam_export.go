package samcmd

import (
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/compgen-io/cgltk/htsio"
	"github.com/spf13/cobra"
)

// columnDef defines a SAM column that can be exported.
type columnDef struct {
	name    string                        // column name (used in header and flag name)
	extract func(rec *htsio.SamRecord) string // extracts the column value
}

var samColumnDefs = []columnDef{
	{"read_name", func(r *htsio.SamRecord) string { return r.ReadName }},
	{"flag", func(r *htsio.SamRecord) string { return fmt.Sprintf("%d", r.Flag) }},
	{"ref_name", func(r *htsio.SamRecord) string { return r.RefName }},
	{"pos", func(r *htsio.SamRecord) string { return fmt.Sprintf("%d", r.Pos) }},
	{"mapq", func(r *htsio.SamRecord) string { return fmt.Sprintf("%d", r.MapQ) }},
	{"cigar", func(r *htsio.SamRecord) string { return r.Cigar }},
	{"ref_next", func(r *htsio.SamRecord) string { return r.RefNext }},
	{"pos_next", func(r *htsio.SamRecord) string { return fmt.Sprintf("%d", r.PosNext) }},
	{"insert_len", func(r *htsio.SamRecord) string { return fmt.Sprintf("%d", r.InsertLen) }},
	{"seq", func(r *htsio.SamRecord) string { return r.Seq }},
	{"qual", func(r *htsio.SamRecord) string { return r.Qual }},
}

var samExportCmd = &cobra.Command{
	GroupID: "samcmd",
	Use:     "sam-export <input.bam>",
	Short:   "Export columns and tags from a SAM/BAM/CRAM file as tab-delimited text",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}

		// Build list of selected columns (in definition order).
		type selectedColumn struct {
			name    string
			extract func(rec *htsio.SamRecord) string
		}
		var selected []selectedColumn

		for i := range samColumnDefs {
			if samExportColumnFlags[i] {
				selected = append(selected, selectedColumn{
					name:    samColumnDefs[i].name,
					extract: samColumnDefs[i].extract,
				})
			}
		}

		tags := parseTags(samExportTags)

		if len(selected) == 0 && len(tags) == 0 {
			return fmt.Errorf("at least one column flag or --tags value is required")
		}

		inputFile := args[0]
		reader, err := htsio.NewSamReader(inputFile)
		if err != nil {
			return err
		}
		defer reader.Close()

		var out io.Writer = os.Stdout
		if samExportOutput != "" && samExportOutput != "-" {
			f, err := os.Create(samExportOutput)
			if err != nil {
				return err
			}
			defer f.Close()
			out = f
		}

		// Write header line.
		headerParts := make([]string, 0, len(selected)+len(tags))
		for _, c := range selected {
			headerParts = append(headerParts, c.name)
		}
		headerParts = append(headerParts, tags...)
		fmt.Fprintln(out, strings.Join(headerParts, "\t"))

		for {
			rec, err := reader.Next()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}

			// Filter by required flags: all bits must be set.
			if samExportRequiredFlag != 0 && rec.Flag&samExportRequiredFlag != samExportRequiredFlag {
				continue
			}
			// Filter by excluded flags: none of these bits may be set.
			if samExportFilterFlag != 0 && rec.Flag&samExportFilterFlag != 0 {
				continue
			}
			// Filter by minimum mapping quality.
			if samExportMinMapQ >= 0 && rec.MapQ < samExportMinMapQ {
				continue
			}

			parts := make([]string, 0, len(selected)+len(tags))
			for _, c := range selected {
				parts = append(parts, c.extract(rec))
			}
			for _, tag := range tags {
				if t, ok := rec.Tags[tag]; ok {
					parts = append(parts, t.Value)
				} else {
					parts = append(parts, "")
				}
			}
			fmt.Fprintln(out, strings.Join(parts, "\t"))
		}

		return nil
	},
}

// parseTags splits the comma-separated tag string into individual tag names.
func parseTags(tagStr string) []string {
	if tagStr == "" {
		return nil
	}
	parts := strings.Split(tagStr, ",")
	var result []string
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p != "" {
			result = append(result, p)
		}
	}
	return result
}

var (
	samExportColumnFlags []bool // one bool per samColumnDefs entry
	samExportTags        string
	samExportOutput      string
	samExportFilterFlag  int
	samExportRequiredFlag int
	samExportMinMapQ     int
)

func init() {
	// Register a --flag for each SAM column.
	samExportColumnFlags = make([]bool, len(samColumnDefs))
	for i := range samColumnDefs {
		name := strings.ReplaceAll(samColumnDefs[i].name, "_", "-")
		samExportCmd.Flags().BoolVar(&samExportColumnFlags[i], name, false, "Include "+samColumnDefs[i].name+" column")
	}

	samExportCmd.Flags().StringVar(&samExportTags, "tags", "", "SAM tags to export (comma-separated, e.g. RX,MI,NM)")
	samExportCmd.Flags().StringVarP(&samExportOutput, "output", "o", "-", "Output file (default: stdout)")
	samExportCmd.Flags().IntVar(&samExportFilterFlag, "filter-flag", 0, "Exclude reads with any of these flag bits set")
	samExportCmd.Flags().IntVar(&samExportRequiredFlag, "required-flag", 0, "Require all of these flag bits to be set")
	samExportCmd.Flags().IntVar(&samExportMinMapQ, "min-mapq", -1, "Minimum mapping quality (exclude reads below this)")
}
