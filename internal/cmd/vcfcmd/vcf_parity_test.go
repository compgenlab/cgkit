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
