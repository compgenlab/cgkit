package seqcmd

import (
	"strings"
	"testing"

	"github.com/compgenlab/cgkit/internal/cmdtest"
)

// Like fastacmd, this package had no tests because its commands printed with
// fmt.Println straight to os.Stdout, where a cobra harness cannot see them.

func run(t *testing.T, args ...string) string {
	t.Helper()
	return cmdtest.Run(t, InitCmd, args...)
}

func TestSeqRevcomp(t *testing.T) {
	cases := []struct{ in, want string }{
		{"ACGT", "ACGT"}, // its own reverse complement
		{"AAAA", "TTTT"},
		{"GATTACA", "TGTAATC"},
		// Case is preserved per base rather than normalized, so a lowercase
		// input stays lowercase and a mixed one stays mixed position by
		// position. Worth pinning: soft-masked reference sequence carries
		// meaning in its case.
		{"gattaca", "tgtaatc"},
		{"aCgT", "AcGt"},
		{"N", "N"},
	}
	for _, c := range cases {
		got := strings.TrimSpace(run(t, "seq-revcomp", c.in))
		if got != c.want {
			t.Errorf("seq-revcomp %s = %q, want %q", c.in, got, c.want)
		}
	}
}

func TestSeqRevcompNeedsExactlyOneArgument(t *testing.T) {
	for _, args := range [][]string{
		{"seq-revcomp"},
		{"seq-revcomp", "ACGT", "TTTT"},
	} {
		if _, err := cmdtest.RunErr(InitCmd, args...); err == nil {
			t.Errorf("%v returned nil, want an argument error", args)
		}
	}
}

// seq-pairwise reports the better of the forward and reverse-complement
// alignments, so an identical pair and a reverse-complemented pair should both
// align cleanly -- the second is the case the reverse pass exists for.
func TestSeqPairwise(t *testing.T) {
	const q = "ACGTACGTAC"

	forward := run(t, "seq-pairwise", q, q)
	if forward == "" {
		t.Fatal("seq-pairwise produced no output")
	}
	if !strings.Contains(forward, q) {
		t.Errorf("an identical pair should align over the whole query:\n%s", forward)
	}

	revcomp := strings.TrimSpace(run(t, "seq-revcomp", q))
	reverse := run(t, "seq-pairwise", q, revcomp)
	if reverse == "" {
		t.Error("seq-pairwise produced no output for a reverse-complemented target")
	}
}

func TestSeqPairwiseNeedsTwoArguments(t *testing.T) {
	for _, args := range [][]string{
		{"seq-pairwise"},
		{"seq-pairwise", "ACGT"},
	} {
		if _, err := cmdtest.RunErr(InitCmd, args...); err == nil {
			t.Errorf("%v returned nil, want an argument error", args)
		}
	}
}

func TestSeqMsa(t *testing.T) {
	got := run(t, "seq-msa", "testdata/simple.fa")
	if got == "" {
		t.Fatal("seq-msa produced no output")
	}
	// Every input sequence has to appear in the alignment, whatever the layout.
	for _, name := range []string{"seq1", "seq2", "seq3"} {
		if !strings.Contains(got, name) {
			t.Errorf("alignment is missing %s:\n%s", name, got)
		}
	}
}

func TestSeqMsaConsensus(t *testing.T) {
	got := run(t, "seq-msa", "--consensus", "testdata/simple.fa")
	if !strings.Contains(got, "consensus") {
		t.Errorf("--consensus produced no consensus row:\n%s", got)
	}
}

// -o must reject a remote locator by name rather than trying to create a local
// file with a scheme in it.
func TestSeqMsaRefusesRemoteOutput(t *testing.T) {
	_, err := cmdtest.RunErr(InitCmd, "seq-msa", "-o", "s3://bucket/out.txt", "testdata/simple.fa")
	if err == nil {
		t.Fatal("a remote -o was accepted")
	}
	if !strings.Contains(err.Error(), "remote locator") {
		t.Errorf("error %q should explain that output must be local", err)
	}
}
