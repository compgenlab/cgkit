package bedcmd

import (
	"strings"
	"testing"
)

// cgkit reads from anywhere and writes locally. Every command taking -o must say
// so by name; these routed through bed.OpenBedWriter with no check, so a remote
// path became an attempt to create a local file literally called
// "s3://bucket/out.bed" and reported a not-found against a path nobody typed.
func TestRemoteOutputIsRefusedByName(t *testing.T) {
	cases := [][]string{
		{"bed-tobed3", "testdata/input.bed", "-o", "s3://bucket/out.bed"},
		{"bed-tobed6", "testdata/input.bed", "-o", "s3://bucket/out.bed"},
		{"bed-clean", "testdata/input.bed", "-o", "s3://bucket/out.bed"},
		{"bed-resize", "-5", "10", "testdata/input.bed", "-o", "s3://bucket/out.bed"},
		{"bed-stats", "testdata/input.bed", "-o", "s3://bucket/out.txt"},
		{"bed-set", "--union", "testdata/setA.bed", "testdata/setB.bed", "-o", "s3://bucket/out.bed"},
	}
	for _, args := range cases {
		_, err := runBedErr(args...)
		if err == nil {
			t.Errorf("%s: a remote -o was accepted", args[0])
			continue
		}
		if !strings.Contains(err.Error(), "remote locator") {
			t.Errorf("%s: error %q should explain that output must be local", args[0], err)
		}
	}
}

// The guard must not swallow stdout, which is what the "apply it after your own
// stdout check" rule existed to prevent. It is now enforced in the function.
func TestStdoutIsNotMistakenForALocator(t *testing.T) {
	for _, out := range []string{"-", ""} {
		if _, err := runBedErr("bed-tobed3", "testdata/input.bed", "-o", out); err != nil {
			t.Errorf("-o %q was rejected: %v", out, err)
		}
	}
}

// The mixed-column-width warning went to os.Stderr, so no test could see it.
// It uses cmd.ErrOrStderr() now, and the harness captures both streams.
func TestBedSetWarnsOnMixedColumnWidths(t *testing.T) {
	// setA is BED6 and bed3.bed is BED3, so strand cannot be honored and the
	// result silently becomes strand-agnostic. Worth saying out loud.
	out, err := runBedErr("bed-set", "--union", "testdata/setA.bed", "testdata/bed3.bed")
	if err != nil {
		t.Fatalf("bed-set: %v", err)
	}
	if !strings.Contains(out, "mixed column widths") {
		t.Errorf("no width warning in output:\n%s", out)
	}
}
