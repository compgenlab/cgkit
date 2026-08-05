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
