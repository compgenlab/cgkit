package bedcmd

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"

	"github.com/spf13/cobra"
)

// resetBedSetFlags restores bed-set's package-global flags to their defaults so
// tests do not leak state into one another (cobra does not reset unset flags).
func resetBedSetFlags() {
	bedSetOutput = "-"
	bedSetInter, bedSetUnion, bedSetSub, bedSetExclusive, bedSetXor = false, false, false, false, false
	bedSetIgnoreStrand, bedSetBed3, bedSetSum, bedSetCount = false, false, false, false
	bedSetDelim, bedSetLabelA, bedSetLabelB = "|", "", ""
	bedSetFlanking = 0
	bedSetTabix = false
}

// runBedErr executes a bed subcommand and returns (stdout, error).
func runBedErr(args ...string) (string, error) {
	root := &cobra.Command{Use: "cgio"}
	InitCmd(root)
	var buf bytes.Buffer
	root.SetOut(&buf)
	root.SetErr(&buf)
	root.SetArgs(args)
	err := root.Execute()
	return buf.String(), err
}

func TestBedSetModes(t *testing.T) {
	const a = "testdata/setA.bed"
	const b = "testdata/setB.bed"
	cases := []struct {
		name string
		args []string
		want string
	}{
		{"inter", []string{"bed-set", "--inter", a, b},
			"chr1\t10\t20\tA1|B1\t0\t+\n"},
		{"union", []string{"bed-set", "--union", a, b},
			"chr1\t0\t30\tA1|B1\t0\t+\n" +
				"chr1\t50\t60\tA2\t0\t+\n" +
				"chr2\t0\t10\tB2\t0\t-\n"},
		{"sub", []string{"bed-set", "--sub", a, b},
			"chr1\t0\t10\tA1\t0\t+\n" +
				"chr1\t50\t60\tA2\t0\t+\n"},
		{"exclusive", []string{"bed-set", "--exclusive", a, b},
			"chr1\t0\t10\tA1\t0\t+\n" +
				"chr1\t20\t30\tB1\t0\t+\n" +
				"chr1\t50\t60\tA2\t0\t+\n" +
				"chr2\t0\t10\tB2\t0\t-\n"},
		{"xor-alias", []string{"bed-set", "--xor", a, b},
			"chr1\t0\t10\tA1\t0\t+\n" +
				"chr1\t20\t30\tB1\t0\t+\n" +
				"chr1\t50\t60\tA2\t0\t+\n" +
				"chr2\t0\t10\tB2\t0\t-\n"},
		{"union-sum", []string{"bed-set", "--union", "--sum", a, b},
			"chr1\t0\t30\tA1|B1\t5\t+\n" +
				"chr1\t50\t60\tA2\t1\t+\n" +
				"chr2\t0\t10\tB2\t5\t-\n"},
		{"union-count-provenance", []string{"bed-set", "--union", "--count", "--a", "setA", "--b", "setB", a, b},
			"chr1\t0\t30\tsetA|setB\t2\t+\n" +
				"chr1\t50\t60\tsetA\t1\t+\n" +
				"chr2\t0\t10\tsetB\t1\t-\n"},
		{"inter-ignore-strand-bed3", []string{"bed-set", "--inter", "--ignore-strand", a, b},
			"chr1\t10\t20\n"},
		{"bed3-flag", []string{"bed-set", "--union", "--bed3", a, b},
			"chr1\t0\t30\n" +
				"chr1\t50\t60\n" +
				"chr2\t0\t10\n"},
		{"union-flanking", []string{"bed-set", "--union", "--flanking", "25", a, b},
			"chr1\t0\t60\tA1|A2|B1\t0\t+\n" +
				"chr2\t0\t10\tB2\t0\t-\n"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			resetBedSetFlags()
			got, err := runBedErr(c.args...)
			if err != nil {
				t.Fatalf("Execute: %v", err)
			}
			if got != c.want {
				t.Errorf("mismatch.\n got: %q\nwant: %q", got, c.want)
			}
		})
	}
}

func TestBedSetModeValidation(t *testing.T) {
	resetBedSetFlags()
	if _, err := runBedErr("bed-set", "testdata/setA.bed", "testdata/setB.bed"); err == nil {
		t.Error("expected error with no mode")
	}
	resetBedSetFlags()
	if _, err := runBedErr("bed-set", "--inter", "--union", "testdata/setA.bed", "testdata/setB.bed"); err == nil {
		t.Error("expected error with two modes")
	}
	resetBedSetFlags()
	if _, err := runBedErr("bed-set", "--inter", "--flanking", "5", "testdata/setA.bed", "testdata/setB.bed"); err == nil {
		t.Error("expected error: --flanking with --inter")
	}
}

func TestBedSetSinks(t *testing.T) {
	dir := t.TempDir()

	// BGZF output (.gz) starts with the gzip/BGZF magic bytes.
	gzPath := filepath.Join(dir, "out.bed.gz")
	resetBedSetFlags()
	if _, err := runBedErr("bed-set", "--union", "-o", gzPath, "testdata/setA.bed", "testdata/setB.bed"); err != nil {
		t.Fatalf("Execute: %v", err)
	}
	raw, err := os.ReadFile(gzPath)
	if err != nil {
		t.Fatal(err)
	}
	if len(raw) < 2 || raw[0] != 0x1f || raw[1] != 0x8b {
		t.Errorf("expected BGZF magic, got % x", raw[:min(2, len(raw))])
	}

	// --tbi writes a companion .tbi index.
	tbiPath := filepath.Join(dir, "idx.bed.gz")
	resetBedSetFlags()
	if _, err := runBedErr("bed-set", "--union", "--tbi", "-o", tbiPath, "testdata/setA.bed", "testdata/setB.bed"); err != nil {
		t.Fatalf("Execute (tbi): %v", err)
	}
	if _, err := os.Stat(tbiPath + ".tbi"); err != nil {
		t.Errorf("expected %s.tbi to exist: %v", tbiPath, err)
	}

	// --tbi to stdout is an error.
	resetBedSetFlags()
	if _, err := runBedErr("bed-set", "--union", "--tbi", "testdata/setA.bed", "testdata/setB.bed"); err == nil {
		t.Error("expected error: --tbi to stdout")
	}
}
