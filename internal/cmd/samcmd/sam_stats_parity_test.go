package samcmd

import (
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/compgenlab/cgkit/internal/cmdtest"
)

// sam-stats is meant as a drop-in replacement for ngsutilsj bam-stats, which
// many pipelines parse, so its report is compared to ngsutilsj's byte for byte.
// The reference binary is found as for the vcf-* parity tests: $NGSUTILSJ, then
// ~/projects/ngsutilsj/dist/ngsutilsj, then ngsutilsj on $PATH. Without one the
// test skips.
func findNgsutilsj(t *testing.T) string {
	t.Helper()
	if p := os.Getenv("NGSUTILSJ"); p != "" {
		if _, err := os.Stat(p); err == nil {
			return p
		}
	}
	home, _ := os.UserHomeDir()
	if cand := filepath.Join(home, "projects", "ngsutilsj", "dist", "ngsutilsj"); fileExists(cand) {
		return cand
	}
	if p, err := exec.LookPath("ngsutilsj"); err == nil {
		return p
	}
	return ""
}

func fileExists(p string) bool { _, err := os.Stat(p); return err == nil }

func TestParityBamStats(t *testing.T) {
	bin := findNgsutilsj(t)
	if bin == "" {
		t.Skip("ngsutilsj reference binary not found; set NGSUTILSJ to enable parity checks")
	}
	gtfPath, samPath := writeGtfFixture(t)
	cases := [][]string{
		{"testdata/test.bam"},
		{"testdata/test_tags.bam"},
		{"--tags", "NH:i,RG:Z,MAPQ", "--show-unmapped-ref", "testdata/test_tags.bam"},
		{"--tags", "NM:i,XS:i,RG:Z,AS:i", "testdata/test_tags.bam"},
		{"--calc-insert", "testdata/test_tags.bam"},
		{"--unique", "testdata/test_tags.bam"},
		{"--rgid", "sample1", "testdata/test_tags.bam"},
		{"--gtf", gtfPath, samPath},
	}
	for _, args := range cases {
		t.Run(args[len(args)-1], func(t *testing.T) {
			// --lenient: these fixtures tag reads with read groups their headers
			// do not declare, which htsjdk otherwise rejects outright. cghts
			// never validates that strictly, and accepts the flag as a no-op.
			java := exec.Command(bin, append([]string{"bam-stats", "--lenient"}, args...)...)
			want, err := java.Output()
			if err != nil {
				t.Fatalf("ngsutilsj bam-stats %v: %v", args, err)
			}
			root, out := cmdtest.NewRoot(InitCmd, append([]string{"sam-stats"}, args...)...)
			root.SetErr(io.Discard) // --gtf announces the GTF load on stderr, as ngsutilsj does
			if err := root.Execute(); err != nil {
				t.Fatalf("sam-stats %v: %v", args, err)
			}
			if got := out.String(); got != string(want) {
				t.Errorf("sam-stats %v differs from ngsutilsj bam-stats:\n--- ngsutilsj\n%s\n--- cgkit\n%s", args, want, got)
			}
		})
	}
}
