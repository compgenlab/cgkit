package vcfcmd

import (
	"os"
	"path/filepath"
	"testing"
)

// brokenVcf is valid through one record and then is not, so a command fails
// partway through the stream with the output file already created.
func brokenVcf(t *testing.T) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "broken.vcf")
	const body = "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tG\t50\tPASS\t.\n" +
		"THIS IS NOT A VCF RECORD\n"
	if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

// A command that dies partway through must not leave its half-written output.
// The writer was never deferred and the close sat after the record loop, so an
// error inside the loop leaked the descriptor and left a file behind -- for a
// bgzip output, one with no BGZF EOF block. That is detectably broken, but only
// if something checks, and meanwhile it sits where a result belongs.
func TestFailedRunLeavesNoPartialOutput(t *testing.T) {
	in := brokenVcf(t)
	cases := []struct {
		name string
		args func(out string) []string
	}{
		{"vcf-clearfilter", func(o string) []string { return []string{"vcf-clearfilter", "-o", o, in} }},
		{"vcf-rename", func(o string) []string {
			return []string{"vcf-rename", "--sample", "A:B", "-o", o, in}
		}},
		{"vcf-chrfix", func(o string) []string { return []string{"vcf-chrfix", "--ucsc", "-o", o, in} }},
		{"vcf-strip", func(o string) []string { return []string{"vcf-strip", "--all", "-o", o, in} }},
		{"vcf-filter", func(o string) []string {
			return []string{"vcf-filter", "--min-qual", "1", "-o", o, in}
		}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			for _, ext := range []string{".vcf", ".vcf.gz"} {
				out := filepath.Join(t.TempDir(), "out"+ext)
				if err := runVcfErr(t, c.args(out)...); err == nil {
					t.Fatalf("%s%s: expected a parse error", c.name, ext)
				}
				if _, err := os.Stat(out); err == nil {
					t.Errorf("%s: a failed run left %s behind", c.name, out)
				}
			}
		})
	}
}

// The counterpart: a run that succeeds keeps its output, so the cleanup above
// cannot be passing by simply never writing anything.
func TestSuccessfulRunKeepsItsOutput(t *testing.T) {
	out := filepath.Join(t.TempDir(), "out.vcf.gz")
	runVcf(t, "vcf-clearfilter", "-o", out, "testdata/sample.vcf")
	st, err := os.Stat(out)
	if err != nil {
		t.Fatalf("a successful run left no output: %v", err)
	}
	if st.Size() == 0 {
		t.Error("output is empty")
	}
}
