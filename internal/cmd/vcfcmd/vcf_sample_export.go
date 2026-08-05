package vcfcmd

import (
	"fmt"
	"io"
	"strings"

	"github.com/compgenlab/cghts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfSampleExportKeys    []string
	vcfSampleExportSamples []string
	vcfSampleExportGT      bool
	vcfSampleExportID      bool
	vcfSampleExportPassing bool
)

var vcfSampleExportCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-sample-export <input.vcf>",
	Short:       "Write sample FORMAT values to a tab-delimited file, one sample per line",
	Long: `Export per-sample FORMAT values as tab-delimited text, one row per sample per
variant. Columns: chrom, pos, [ID], ref, alt, sample, then each exported key.

  --key, -k KEY      FORMAT key to export (glob, repeatable, at least one required)
  --sample, -s NAME  sample to export (glob, repeatable; default all)
  --gt               export GT and convert "0/1" to ref/alt bases
  --id               include the ID column
  --passing          only export passing variants`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		validKeys := append([]string(nil), vcfSampleExportKeys...)
		convertGT := vcfSampleExportGT
		if convertGT {
			validKeys = append(validKeys, "GT")
		}
		if len(validKeys) == 0 {
			return fmt.Errorf("you must specify at least one field to export (--key or --gt)")
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

		samples := selectByGlob(header.Samples(), vcfSampleExportSamples)
		keys := selectByGlob(header.FormatIDs(), validKeys)

		// selectByGlob filters against the header's declared FORMAT ids, so a
		// --key naming a field the header never declared -- common in
		// hand-assembled or tool-truncated VCFs -- yielded no column and no
		// word about it. The "at least one field" guard above has already
		// passed by then, so the result was a table of locus columns and
		// nothing else.
		warnUnmatchedKeys(cmd, validKeys, keys)

		out := cmd.OutOrStdout()
		var hdr []string
		hdr = append(hdr, "chrom", "pos")
		if vcfSampleExportID {
			hdr = append(hdr, "ID")
		}
		hdr = append(hdr, "ref", "alt", "sample")
		hdr = append(hdr, keys...)
		fmt.Fprintln(out, strings.Join(hdr, "\t"))

		for {
			rec, err := reader.NextRecord()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if vcfSampleExportPassing && rec.IsFiltered() {
				continue
			}
			for _, s := range samples {
				idx := header.SampleIndex(s)
				row := []string{rec.Chrom, fmt.Sprint(rec.Pos)}
				if vcfSampleExportID {
					row = append(row, idOrMissing(rec))
				}
				row = append(row, rec.Ref, rec.AltOrig(), s)
				attr, aerr := rec.Sample(idx)
				if aerr != nil {
					return aerr
				}
				for _, k := range keys {
					row = append(row, sampleExportValue(rec, attr, idx, k, convertGT))
				}
				fmt.Fprintln(out, strings.Join(row, "\t"))
			}
		}
		return nil
	},
}

// selectByGlob returns the candidates that exactly match or glob-match any of
// the patterns, preserving candidate order and de-duplicating. A nil/empty
// pattern list selects every candidate.
func selectByGlob(candidates, patterns []string) []string {
	if len(patterns) == 0 {
		return append([]string(nil), candidates...)
	}
	var out []string
	seen := map[string]bool{}
	for _, c := range candidates {
		match := false
		for _, p := range patterns {
			if c == p || globMatch(c, p) {
				match = true
				break
			}
		}
		if match && !seen[c] {
			seen[c] = true
			out = append(out, c)
		}
	}
	return out
}

func idOrMissing(rec *vcf.VcfRecord) string {
	if id := rec.ID(); id != "" {
		return id
	}
	return "."
}

func sampleExportValue(rec *vcf.VcfRecord, attr *vcf.Attributes, idx int, key string, convertGT bool) string {
	if key == "GT" && convertGT {
		if gt, ok := rec.GenotypeBases(idx); ok {
			return gt
		}
		return "."
	}
	if v, ok := attr.Get(key); ok {
		return v.String()
	}
	return "."
}

func init() {
	f := vcfSampleExportCmd.Flags()
	f.StringArrayVarP(&vcfSampleExportKeys, "key", "k", nil, "FORMAT key to export (glob, repeatable)")
	f.StringArrayVarP(&vcfSampleExportSamples, "sample", "s", nil, "Sample to export (glob, repeatable; default all)")
	f.BoolVar(&vcfSampleExportGT, "gt", false, "Export GT and convert to ref/alt bases")
	f.BoolVar(&vcfSampleExportID, "id", false, "Include the ID column")
	f.BoolVar(&vcfSampleExportPassing, "passing", false, "Only export passing variants")
}

// warnUnmatchedKeys reports requested keys that matched no declared FORMAT id.
//
// A warning rather than an error: a glob that matches nothing is a reasonable
// thing to pass across a set of files with differing headers, and the export is
// still valid for the keys that did match. Silence was the problem.
func warnUnmatchedKeys(cmd *cobra.Command, requested, matched []string) {
	if len(requested) == 0 {
		return
	}
	got := make(map[string]bool, len(matched))
	for _, k := range matched {
		got[k] = true
	}
	var missing []string
	for _, want := range requested {
		if got[want] {
			continue
		}
		// A glob is satisfied by anything it matched, so only report one that
		// matched nothing at all.
		if isGlobPattern(want) && len(matched) > 0 && globMatchedAny(want, matched) {
			continue
		}
		missing = append(missing, want)
	}
	if len(missing) == 0 {
		return
	}
	fmt.Fprintf(cmd.ErrOrStderr(),
		"warning: no FORMAT field declared in the header matches: %s\n",
		strings.Join(missing, ", "))
	if len(matched) == 0 {
		fmt.Fprintln(cmd.ErrOrStderr(),
			"         the output will carry no per-sample columns")
	}
}

// isGlobPattern reports whether a requested key is a pattern rather than a
// literal id.
func isGlobPattern(s string) bool {
	return strings.ContainsAny(s, "*?[")
}

// globMatchedAny reports whether a pattern accounts for at least one match.
func globMatchedAny(pattern string, matched []string) bool {
	for _, m := range matched {
		if globMatch(m, pattern) {
			return true
		}
	}
	return false
}
