package vcfcmd

import (
	"testing"
)

// TestParityGroup3 checks the SV converters vcf-tobedpe and vcf-svtofasta
// against the ngsutilsj reference. Both produce non-VCF output, so output is
// compared directly after dropping "##" comment lines.
func TestParityGroup3(t *testing.T) {
	bin := findNgsutilsj(t)
	if bin == "" {
		t.Skip("ngsutilsj reference binary not found; set NGSUTILSJ to enable parity checks")
	}
	const sv = "testdata/sv.vcf"
	const ref = "testdata/svref.fa"

	// vcf-tobedpe variants (cgio uses the same flag names as ngsutilsj).
	bedpeCases := []struct {
		name string
		args []string
	}{
		{"plain", []string{"vcf-tobedpe", sv}},
		{"no-del-offset", []string{"vcf-tobedpe", "--no-del-offset", sv}},
		{"name-id", []string{"vcf-tobedpe", "--name", "@ID", sv}},
		{"name-info", []string{"vcf-tobedpe", "--name", "SVTYPE", sv}},
		// A bare --score key uses the INFO field. (cgio also accepts an explicit
		// ":INFO", which ngsutilsj rejects — so parity uses the bare form.)
		{"score-info", []string{"vcf-tobedpe", "--score", "SVTYPE", sv}},
		{"unique-event", []string{"vcf-tobedpe", "--unique-event", sv}},
	}
	for _, tc := range bedpeCases {
		t.Run("tobedpe-"+tc.name, func(t *testing.T) {
			want := stripComments(runJava(t, bin, tc.args...))
			got := stripComments(runVcf(t, tc.args...))
			if got != want {
				t.Errorf("tobedpe %s parity\n java: %q\n cgio: %q", tc.name, want, got)
			}
		})
	}

	// vcf-svtofasta: cgio's --bnd corresponds to ngsutilsj's (typo'd) --bnf.
	svCases := []struct {
		name      string
		cgioExtra []string
		javaExtra []string
	}{
		{"bnd", []string{"--bnd"}, []string{"--bnf"}},
		{"include-ref", []string{"--bnd", "--include-ref"}, []string{"--bnf", "--include-ref"}},
	}
	for _, tc := range svCases {
		t.Run("svtofasta-"+tc.name, func(t *testing.T) {
			cgioArgs := append([]string{"vcf-svtofasta"}, tc.cgioExtra...)
			cgioArgs = append(cgioArgs, "--flanking", "10", ref, sv)
			javaArgs := append([]string{"vcf-svtofasta"}, tc.javaExtra...)
			javaArgs = append(javaArgs, "--flanking", "10", ref, sv)
			want := runJava(t, bin, javaArgs...)
			got := runVcf(t, cgioArgs...)
			if got != want {
				t.Errorf("svtofasta %s parity\n java: %q\n cgio: %q", tc.name, want, got)
			}
		})
	}
}
