package samcmd

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cgkit/internal/cmdtest"
)

// sam-tofasta and sam-tofastq had no tests, because both wrote to os.Stdout
// where a cobra harness cannot see them. They use cmd.OutOrStdout() now.

// tofastxSAM has one forward read and one reverse-strand read. The second is
// the interesting one: a reverse-strand alignment stores the sequence against
// the reference, so emitting it verbatim would give the wrong strand of the
// original read.
const tofastxSAM = "@HD\tVN:1.6\n" +
	"@SQ\tSN:chr1\tLN:1000\n" +
	"fwd\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tBC:Z:AAAA\n" +
	"rev\t16\tchr1\t20\t60\t4M\t*\t0\t0\tAAAC\tJJJJ\tBC:Z:CCCC\n"

func tofastxInput(t *testing.T) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "in.sam")
	if err := os.WriteFile(p, []byte(tofastxSAM), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestSamToFasta(t *testing.T) {
	got := cmdtest.Run(t, InitCmd, "sam-tofasta", tofastxInput(t))
	want := ">fwd\nACGT\n" +
		// AAAC reverse-complemented back to the read's own orientation.
		">rev\nGTTT\n"
	if got != want {
		t.Errorf("sam-tofasta mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestSamToFastq(t *testing.T) {
	got := cmdtest.Run(t, InitCmd, "sam-tofastq", tofastxInput(t))
	if !strings.Contains(got, "@fwd\nACGT\n+\nIIII\n") {
		t.Errorf("forward read wrong:\n%s", got)
	}
	// The quality string is reversed alongside the sequence, or the two no
	// longer line up base for base.
	if !strings.Contains(got, "@rev\nGTTT\n+\nJJJJ\n") {
		t.Errorf("reverse read wrong:\n%s", got)
	}
}

// --write-tag puts SAM tags in the record comment, which is how read-level metadata
// survives the conversion at all.
func TestSamToFastqTags(t *testing.T) {
	got := cmdtest.Run(t, InitCmd, "sam-tofastq", "--write-tag", "BC", tofastxInput(t))
	if !strings.Contains(got, "BC:Z:AAAA") || !strings.Contains(got, "BC:Z:CCCC") {
		t.Errorf("--write-tag BC did not reach the output:\n%s", got)
	}
}

func TestSamToFastxRefusesRemoteOutput(t *testing.T) {
	in := tofastxInput(t)
	for _, name := range []string{"sam-tofasta", "sam-tofastq"} {
		_, err := cmdtest.RunErr(InitCmd, name, in, "s3://bucket/out.fa")
		if err == nil {
			t.Errorf("%s accepted a remote output", name)
			continue
		}
		if !strings.Contains(err.Error(), "remote locator") {
			t.Errorf("%s: error %q should explain that output must be local", name, err)
		}
	}
}

func TestSamToFastxNeedsAnInput(t *testing.T) {
	for _, name := range []string{"sam-tofasta", "sam-tofastq"} {
		if _, err := cmdtest.RunErr(InitCmd, name); err == nil {
			t.Errorf("%s with no arguments returned nil", name)
		}
	}
}
