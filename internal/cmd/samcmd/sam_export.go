package samcmd

import (
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/compgen-io/cgkit/htsio"
	_ "github.com/compgen-io/cgkit/htsio/bam"
	_ "github.com/compgen-io/cgkit/htsio/cram"
	_ "github.com/compgen-io/cgkit/htsio/sam"
	"github.com/spf13/cobra"
)

// columnDef defines a SAM column that can be exported.
type columnDef struct {
	name    string                             // column name (used in header and flag name)
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

		opts, err := samExportReaderFlags.buildReaderOpts()
		if err != nil {
			return err
		}

		inputFile := args[0]
		reader, err := htsio.NewSamReader(inputFile, opts)
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

		for rec, err := range reader.Records() {
			if err != nil {
				return err
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
	samExportColumnFlags  []bool // one bool per samColumnDefs entry
	samExportTags         string
	samExportOutput       string
	samExportReaderFlags  samReaderFlags
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

	samExportReaderFlags.register(samExportCmd)
}
