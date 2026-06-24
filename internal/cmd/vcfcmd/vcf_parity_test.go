package vcfcmd

import (
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// findNgsutilsj locates the ngsutilsj reference binary, used to verify output
// parity. It returns "" (and the caller skips) when the binary is not present,
// so this test is a no-op in environments without the reference (e.g. CI).
func findNgsutilsj(t *testing.T) string {
	t.Helper()
	if p := os.Getenv("NGSUTILSJ"); p != "" {
		if _, err := os.Stat(p); err == nil {
			return p
		}
	}
	home, _ := os.UserHomeDir()
	cand := filepath.Join(home, "projects", "ngsutilsj", "dist", "ngsutilsj")
	if _, err := os.Stat(cand); err == nil {
		return cand
	}
	if p, err := exec.LookPath("ngsutilsj"); err == nil {
		return p
	}
	return ""
}

func runJava(t *testing.T, bin string, args ...string) string {
	t.Helper()
	cmd := exec.Command(bin, args...)
	out, err := cmd.Output() // stderr (progress) is discarded
	if err != nil {
		t.Fatalf("ngsutilsj %v: %v", args, err)
	}
	return string(out)
}

// stripProvenance drops the non-deterministic ##fileDate and ##ngsutilsj_*
// provenance header lines that cgio intentionally does not emit.
func stripProvenance(s string) string {
	var keep []string
	for _, line := range strings.Split(s, "\n") {
		if strings.HasPrefix(line, "##fileDate") || strings.HasPrefix(line, "##ngsutilsj") || strings.HasPrefix(line, "##cgio") {
			continue
		}
		keep = append(keep, line)
	}
	return strings.Join(keep, "\n")
}

// TestParityWithNgsutilsj verifies that cgio's output matches the ngsutilsj
// reference for the commands that are designed to be byte-identical (after
// removing provenance header lines cgio omits by design).
func TestParityWithNgsutilsj(t *testing.T) {
	bin := findNgsutilsj(t)
	if bin == "" {
		t.Skip("ngsutilsj reference binary not found; set NGSUTILSJ to enable parity checks")
	}
	const vcf = "testdata/sample.vcf"

	cases := []struct {
		name string
		args []string
	}{
		{"samples", []string{"vcf-samples", vcf}},
		{"tobed", []string{"vcf-tobed", vcf}},
		{"tobed-passing", []string{"vcf-tobed", "--passing", vcf}},
		{"tobed-includepos-pad", []string{"vcf-tobed", "--include-pos", "--padding", "5", vcf}},
		{"stats", []string{"vcf-stats", vcf}},
		{"stats-info", []string{"vcf-stats", "--info-tally", "SVTYPE", "--info-present", "DB", vcf}},
		{"stats-filtercombo", []string{"vcf-stats", "--filter-combo", vcf}},
		{"tstv", []string{"vcf-tstv", vcf}},
		{"tstv-passing", []string{"vcf-tstv", "--passing", vcf}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want := stripProvenance(runJava(t, bin, tc.args...))
			got := stripProvenance(runVcf(t, tc.args...))
			if got != want {
				t.Errorf("parity mismatch for %v\n java: %q\n cgio: %q", tc.args, want, got)
			}
		})
	}
}

// dataRows returns only the non-header (non-#) lines of VCF/tab output.
func dataRows(s string) string {
	var keep []string
	for _, line := range strings.Split(strings.TrimRight(s, "\n"), "\n") {
		if !strings.HasPrefix(line, "#") {
			keep = append(keep, line)
		}
	}
	return strings.Join(keep, "\n")
}

// TestParityAnnotate verifies vcf-annotate output against ngsutilsj. Group A
// annotators, copy-logratio, and vardist are byte-identical on the data rows.
func TestParityAnnotate(t *testing.T) {
	bin := findNgsutilsj(t)
	if bin == "" {
		t.Skip("ngsutilsj reference binary not found; set NGSUTILSJ to enable parity checks")
	}
	const vcf = "testdata/annotate.vcf"
	cases := []struct {
		name string
		args []string
	}{
		{"groupA", []string{"vcf-annotate", "--indel", "--tstv", "--dosage", "--auto-id", vcf}},
		{"copy-logratio", []string{"vcf-annotate", "--copy-logratio", "TUMOR:NORMAL", vcf}},
		{"vardist", []string{"vcf-annotate", "--vardist", vcf}},
		{"single-tag", []string{"vcf-annotate", "--tag", "PANEL:myset", vcf}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want := dataRows(runJava(t, bin, tc.args...))
			got := dataRows(runVcf(t, tc.args...))
			if got != want {
				t.Errorf("annotate parity (%s)\n java: %q\n cgio: %q", tc.name, want, got)
			}
		})
	}
}

// normQual strips a trailing ".0" from the QUAL column (field 5) of each data
// row. cgio writes records it did not modify verbatim (preserving "20.0"),
// whereas ngsutilsj rebuilds every record ("20"); this isolates that
// intentional difference so the annotation values can be compared.
func normQual(s string) string {
	var out []string
	for _, line := range strings.Split(dataRows(s), "\n") {
		f := strings.Split(line, "\t")
		if len(f) > 5 {
			f[5] = strings.TrimSuffix(f[5], ".0")
		}
		out = append(out, strings.Join(f, "\t"))
	}
	return strings.Join(out, "\n")
}

// TestParityAnnotateBedTab verifies the BED/tabix annotators against ngsutilsj.
// Both read the same bgzipped+indexed files. Annotation values match; QUAL on
// untouched rows is normalized (see normQual).
func TestParityAnnotateBedTab(t *testing.T) {
	bin := findNgsutilsj(t)
	if bin == "" {
		t.Skip("ngsutilsj reference binary not found; set NGSUTILSJ to enable parity checks")
	}
	const vcf = "testdata/annotate.vcf"
	cases := []struct {
		name string
		args []string
	}{
		{"bed", []string{"vcf-annotate", "--bed", "REGION:testdata/regions.bed.gz", vcf}},
		{"bed-flag", []string{"vcf-annotate", "--bed-flag", "INREG:testdata/regions.bed.gz", vcf}},
		{"tab-num", []string{"vcf-annotate", "--tab", "SCORE:testdata/scores.tab.gz,5,n", vcf}},
		{"tab-alt", []string{"vcf-annotate", "--tab", "LBL:testdata/scores.tab.gz,6,alt=4", vcf}},
		{"tab-colname", []string{"vcf-annotate", "--tab", "LBL:testdata/scores_hdr.tab.gz,label,alt=alt,ref=ref", vcf}},
		{"vcf-field", []string{"vcf-annotate", "--vcf", "KAF:AF:testdata/source.vcf.gz", vcf}},
		{"vcf-passing", []string{"vcf-annotate", "--vcf", "KAF:AF:testdata/source.vcf.gz:@", vcf}},
		{"vcf-exact", []string{"vcf-annotate", "--vcf", "KAF:AF:testdata/source.vcf.gz:!", vcf}},
		{"vcf-id", []string{"vcf-annotate", "--vcf-id", "testdata/source.vcf.gz", vcf}},
		{"vcf-flag", []string{"vcf-annotate", "--vcf-flag", "KNOWN:testdata/source.vcf.gz", vcf}},
		{"in-file-flag", []string{"vcf-annotate", "--in-file", "HIT:DP:testdata/dp_set.txt", vcf}},
		{"in-file-tabcol", []string{"vcf-annotate", "--in-file", "LABEL:DP:testdata/dp_val.txt:tabcol=2", vcf}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want := normQual(runJava(t, bin, tc.args...))
			got := normQual(runVcf(t, tc.args...))
			if got != want {
				t.Errorf("bed/tab parity (%s)\n java: %q\n cgio: %q", tc.name, want, got)
			}
		})
	}
}

// TestParityAnnotateFlanking verifies --flanking against ngsutilsj on a small
// reference. Both read the same indexed ref.fa.
func TestParityAnnotateFlanking(t *testing.T) {
	bin := findNgsutilsj(t)
	if bin == "" {
		t.Skip("ngsutilsj reference binary not found; set NGSUTILSJ to enable parity checks")
	}
	args := []string{"vcf-annotate", "--flanking", "testdata/ref.fa", "testdata/flank.vcf"}
	want := dataRows(runJava(t, bin, args...))
	got := dataRows(runVcf(t, args...))
	if got != want {
		t.Errorf("flanking parity\n java: %q\n cgio: %q", want, got)
	}
}

// TestParityAnnotateGroupBValues verifies the sample-count annotators produce
// the same values; the FORMAT column ordering differs by design (cgio uses a
// stable order), so the per-line tokens are compared order-insensitive.
func TestParityAnnotateGroupBValues(t *testing.T) {
	bin := findNgsutilsj(t)
	if bin == "" {
		t.Skip("ngsutilsj reference binary not found; set NGSUTILSJ to enable parity checks")
	}
	args := []string{"vcf-annotate", "--vaf", "--minor-strand", "--fisher-sb", "testdata/annotate.vcf"}
	normalize := func(s string) string {
		var lines []string
		for _, line := range strings.Split(dataRows(s), "\n") {
			toks := strings.FieldsFunc(line, func(r rune) bool { return r == '\t' || r == ':' })
			sort.Strings(toks)
			lines = append(lines, strings.Join(toks, " "))
		}
		return strings.Join(lines, "\n")
	}
	want := normalize(runJava(t, bin, args...))
	got := normalize(runVcf(t, args...))
	if got != want {
		t.Errorf("annotate group-B value parity\n java: %q\n cgio: %q", want, got)
	}
}

// TestParityFilter verifies vcf-filter against ngsutilsj. The FILTER codes are a
// set, so their order is normalized (cgio applies filters in command-line order;
// ngsutilsj uses its option-invocation order), along with QUAL on untouched rows.
func TestParityFilter(t *testing.T) {
	bin := findNgsutilsj(t)
	if bin == "" {
		t.Skip("ngsutilsj reference binary not found; set NGSUTILSJ to enable parity checks")
	}
	norm := func(s string) string {
		var lines []string
		for _, line := range strings.Split(dataRows(s), "\n") {
			f := strings.Split(line, "\t")
			if len(f) > 5 {
				f[5] = strings.TrimSuffix(f[5], ".0") // QUAL
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
	const vcf = "testdata/annotate.vcf"
	cases := []struct {
		name string
		args []string
	}{
		{"snv-qual-chrom", []string{"vcf-filter", "--snv", "--qual", "30", "--chrom-fail", "chr2", vcf}},
		{"indel-passing", []string{"vcf-filter", "--indel", "--passing", vcf}},
		{"het-hom", []string{"vcf-filter", "--het", "--hom", vcf}},
		{"maxdel-failing", []string{"vcf-filter", "--max-del", "1", "--failing", vcf}},
		{"chrom-pass", []string{"vcf-filter", "--chrom-pass", "chr1", vcf}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want := norm(runJava(t, bin, tc.args...))
			got := norm(runVcf(t, tc.args...))
			if got != want {
				t.Errorf("filter parity (%s)\n java: %q\n cgio: %q", tc.name, want, got)
			}
		})
	}
}

// TestParityExportValues verifies that cgio and ngsutilsj produce the same set
// of exported values. The column ordering differs by design (cgio uses a stable
// order), so the comparison is order-insensitive per line.
func TestParityExportValues(t *testing.T) {
	bin := findNgsutilsj(t)
	if bin == "" {
		t.Skip("ngsutilsj reference binary not found; set NGSUTILSJ to enable parity checks")
	}
	args := []string{"vcf-export", "--no-vcf-header", "--id", "--qual", "--filter",
		"--info", "DP", "--info", "AF", "--format", "AD", "testdata/sample.vcf"}
	normalize := func(s string) string {
		var lines []string
		for _, line := range strings.Split(strings.TrimRight(s, "\n"), "\n") {
			toks := strings.Split(line, "\t")
			sort.Strings(toks)
			lines = append(lines, strings.Join(toks, "\t"))
		}
		return strings.Join(lines, "\n")
	}
	want := normalize(runJava(t, bin, args...))
	got := normalize(runVcf(t, args...))
	if got != want {
		t.Errorf("export value parity mismatch\n java: %q\n cgio: %q", want, got)
	}
}
