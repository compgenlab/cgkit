package fastacmd

import (
	"strings"
	"testing"

	"github.com/compgenlab/cgkit/internal/cmdtest"
)

// This package had no tests at all, and the reason was mechanical rather than
// deliberate: every command wrote to os.Stdout, which a cobra harness cannot
// capture. They write to cmd.OutOrStdout() now, so here is the first coverage
// they have had.

func run(t *testing.T, args ...string) string {
	t.Helper()
	return cmdtest.Run(t, InitCmd, args...)
}

func TestFastaCat(t *testing.T) {
	got := run(t, "fasta-cat", "testdata/simple.fa")
	// Unwrapped, one line per sequence, comments preserved on the header.
	want := ">seq1 first record\n" +
		"ACGTACGTACGGCCAATTGG\n" +
		">seq2\n" +
		"AAAATTTT\n" +
		">seq3 gc-rich\n" +
		"GGGGCCCC\n"
	if got != want {
		t.Errorf("fasta-cat mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestFastaWrap(t *testing.T) {
	got := run(t, "fasta-wrap", "-w", "5", "testdata/simple.fa")
	// seq1 is 20 bases, so at width 5 it is exactly four lines.
	if !strings.Contains(got, "ACGTA\nCGTAC\nGGCCA\nATTGG\n") {
		t.Errorf("sequence was not wrapped at 5:\n%s", got)
	}
	for _, line := range strings.Split(strings.TrimRight(got, "\n"), "\n") {
		if strings.HasPrefix(line, ">") {
			continue
		}
		if len(line) > 5 {
			t.Errorf("line %q exceeds the requested width", line)
		}
	}
}

// -1 means "do not wrap", which is a different thing from a width of zero and
// is easy to get wrong when the flag is clamped.
func TestFastaWrapNoWrapping(t *testing.T) {
	got := run(t, "fasta-wrap", "-w", "-1", "testdata/simple.fa")
	if !strings.Contains(got, "ACGTACGTACGGCCAATTGG\n") {
		t.Errorf("-w -1 should emit the sequence on one line:\n%s", got)
	}
}

func TestFastaGC(t *testing.T) {
	got := run(t, "fasta-gc", "testdata/simple.fa")
	want := "seq1\t0.5500\n" + // 11 G/C of 20
		"seq2\t0.0000\n" + // no G or C
		"seq3\t1.0000\n" // all G/C
	if got != want {
		t.Errorf("fasta-gc mismatch.\n got: %q\nwant: %q", got, want)
	}
}

// Every command in this group takes an input file, so none of them should treat
// a missing argument as a reason to succeed.
func TestNoArgsIsAnError(t *testing.T) {
	for _, name := range []string{"fasta-cat", "fasta-wrap", "fasta-gc"} {
		if _, err := cmdtest.RunErr(InitCmd, name); err == nil {
			t.Errorf("%s with no arguments returned nil", name)
		}
	}
}

// InitCmd must not depend on fastqcmd having run first, or either package alone
// registers commands against a group cobra does not know about.
func TestInitCmdStandsAlone(t *testing.T) {
	root, _ := cmdtest.NewRoot(InitCmd)
	if !root.ContainsGroup("fastaqcmd") {
		t.Error("fastacmd.InitCmd did not register its own help group")
	}
}
