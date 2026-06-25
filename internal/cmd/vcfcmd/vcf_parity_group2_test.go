package vcfcmd

import (
	"path/filepath"
	"testing"
)

// TestParityGroup2 checks vcf-concat (default and --chunks) and vcf-merge
// against the ngsutilsj reference, comparing normalized data rows.
func TestParityGroup2(t *testing.T) {
	bin := findNgsutilsj(t)
	if bin == "" {
		t.Skip("ngsutilsj reference binary not found; set NGSUTILSJ to enable parity checks")
	}

	t.Run("concat", func(t *testing.T) {
		args := []string{"vcf-concat", "testdata/concat_a.vcf", "testdata/concat_b.vcf"}
		want := normVcfData(runJava(t, bin, args...))
		got := normVcfData(runVcf(t, args...))
		if got != want {
			t.Errorf("concat parity\n java: %q\n cgio: %q", want, got)
		}
	})

	t.Run("merge", func(t *testing.T) {
		args := []string{"vcf-merge", "testdata/merge_a.vcf", "testdata/merge_b.vcf"}
		want := normVcfData(runJava(t, bin, args...))
		got := normVcfData(runVcf(t, args...))
		if got != want {
			t.Errorf("merge parity\n java: %q\n cgio: %q", want, got)
		}
	})

	t.Run("concat-chunks", func(t *testing.T) {
		// Split sample.vcf into numbered chunks; cgio --chunks must match
		// ngsutilsj vcf-concat-n on the same .1 starting file.
		base := filepath.Join(t.TempDir(), "split")
		runVcf(t, "vcf-split", "--out", base, "--num", "2", "testdata/sample.vcf")
		first := base + ".1.vcf.gz"

		want := normVcfData(runJava(t, bin, "vcf-concat-n", first))
		got := normVcfData(runVcf(t, "vcf-concat", "--chunks", first))
		if got != want {
			t.Errorf("concat --chunks parity\n java: %q\n cgio: %q", want, got)
		}
	})
}
