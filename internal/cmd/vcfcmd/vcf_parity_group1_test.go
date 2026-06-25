package vcfcmd

import (
	"sort"
	"strings"
	"testing"
)

// normVcfData reduces VCF output to its data rows with QUAL ".0" trimmed and
// FILTER codes sorted, matching the documented deviations (cgio keeps untouched
// QUALs verbatim and emits FILTER codes in command-line order).
func normVcfData(s string) string {
	var lines []string
	for _, line := range strings.Split(dataRows(s), "\n") {
		f := strings.Split(line, "\t")
		if len(f) > 5 {
			f[5] = strings.TrimSuffix(f[5], ".0")
		}
		if len(f) > 6 {
			codes := strings.Split(f[6], ";")
			sort.Strings(codes)
			f[6] = strings.Join(codes, ";")
		}
		lines = append(lines, strings.Join(f, "\t"))
	}
	return strings.Join(lines, "\n")
}

// stripComments drops blank and "##"-prefixed lines (tab outputs whose only
// differing lines are cgio-branded provenance comments).
func stripComments(s string) string {
	var keep []string
	for _, line := range strings.Split(strings.TrimRight(s, "\n"), "\n") {
		if line == "" || strings.HasPrefix(line, "##") {
			continue
		}
		keep = append(keep, line)
	}
	return strings.Join(keep, "\n")
}

func chromLine(s string) string {
	for _, line := range strings.Split(s, "\n") {
		if strings.HasPrefix(line, "#CHROM") {
			return line
		}
	}
	return ""
}

// TestParityGroup1 checks the group-1 commands against the ngsutilsj reference.
func TestParityGroup1(t *testing.T) {
	bin := findNgsutilsj(t)
	if bin == "" {
		t.Skip("ngsutilsj reference binary not found; set NGSUTILSJ to enable parity checks")
	}
	const sample = "testdata/sample.vcf"
	const annotate = "testdata/annotate.vcf"

	// VCF-output commands: compare normalized data rows.
	vcfCases := []struct {
		name string
		args []string
	}{
		{"clearfilter", []string{"vcf-clearfilter", "--filter", "lowqual", sample}},
		{"chrfix-ensembl", []string{"vcf-chrfix", "--ensembl", sample}},
		{"chrfix-contigs", []string{"vcf-chrfix", "--ensembl", "--contigs", "1", sample}},
		{"remove-flags", []string{"vcf-remove-flags", sample}},
		{"strip-all", []string{"vcf-strip", "--all", sample}},
		{"strip-keep-info", []string{"vcf-strip", "--info", "*", "--keep-info", "DP", sample}},
		{"strip-dbsnp", []string{"vcf-strip", "--dbsnp", sample}},
		{"strip-format", []string{"vcf-strip", "--format", "AD", sample}},
		{"strip-only-snvs", []string{"vcf-strip", "--info", "*", "--only-snvs", sample}},
	}
	for _, tc := range vcfCases {
		t.Run(tc.name, func(t *testing.T) {
			want := normVcfData(runJava(t, bin, tc.args...))
			got := normVcfData(runVcf(t, tc.args...))
			if got != want {
				t.Errorf("%s parity\n java: %q\n cgio: %q", tc.name, want, got)
			}
		})
	}

	// vcf-rename only changes the #CHROM sample columns.
	t.Run("rename", func(t *testing.T) {
		// ngsutilsj's renameSample only supports renaming by name (a numeric
		// old-name throws), so parity uses name-based renames; cgio additionally
		// accepts a sample number.
		args := []string{"vcf-rename", "--sample", "NORMAL:GERMLINE", "--sample", "TUMOR:CANCER", sample}
		if got, want := chromLine(runVcf(t, args...)), chromLine(runJava(t, bin, args...)); got != want {
			t.Errorf("rename #CHROM parity\n java: %q\n cgio: %q", want, got)
		}
	})

	// Tab-output commands: compare non-comment lines directly.
	tabCases := []struct {
		name string
		args []string
	}{
		{"header-info-info", []string{"vcf-header-info", "--info", sample}},
		{"header-info-format", []string{"vcf-header-info", "--format", sample}},
		{"header-info-filters", []string{"vcf-header-info", "--filters", sample}},
		{"header-info-contig", []string{"vcf-header-info", "--contig", sample}},
		{"sample-export", []string{"vcf-sample-export", "--gt", annotate}},
		{"tocount", []string{"vcf-tocount", "--af", "--total", annotate}},
		{"tocount-roao", []string{"vcf-tocount", "--use-ro-ao", "--af", "testdata/roao.vcf"}},
	}
	for _, tc := range tabCases {
		t.Run(tc.name, func(t *testing.T) {
			want := stripComments(runJava(t, bin, tc.args...))
			got := stripComments(runVcf(t, tc.args...))
			if got != want {
				t.Errorf("%s parity\n java: %q\n cgio: %q", tc.name, want, got)
			}
		})
	}
}
