package bedcmd

import (
	"bytes"
	"testing"

	"github.com/spf13/cobra"
)

// runBed executes a bed subcommand against a fresh root and returns its stdout.
func runBed(t *testing.T, args ...string) string {
	t.Helper()
	root := &cobra.Command{Use: "cgkit"}
	InitCmd(root)
	var buf bytes.Buffer
	root.SetOut(&buf)
	root.SetErr(&buf)
	root.SetArgs(args)
	if err := root.Execute(); err != nil {
		t.Fatalf("Execute(%v): %v", args, err)
	}
	return buf.String()
}

func TestBedToBed3(t *testing.T) {
	want := "chr1\t100\t200\n" +
		"chr2\t50\t75\n" +
		"chr3\t10\t20\n" +
		"chr1\t100\t200\n"
	if got := runBed(t, "bed-tobed3", "testdata/input.bed"); got != want {
		t.Errorf("bed-tobed3 mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestBedToBed6(t *testing.T) {
	// Strict BED6: extras dropped, BED3 line expanded, float scores retained.
	want := "chr1\t100\t200\t\t0\t+\n" +
		"chr2\t50\t75\tregionB\t2.9\t-\n" +
		"chr3\t10\t20\tregionC\t0\t+\n" +
		"chr1\t100\t200\tregionA\t1.9\t+\n"
	if got := runBed(t, "bed-tobed6", "testdata/input.bed"); got != want {
		t.Errorf("bed-tobed6 mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestBedResize(t *testing.T) {
	// -5 100: plus shifts start left, minus extends end, none clamps start at 0
	// (strand "." is forced to "+" on output).
	want := "chr1\t900\t2000\tplus\t5\t+\n" +
		"chr1\t3000\t4100\tminus\t10\t-\n" +
		"chr2\t0\t150\tnone\t0\t+\n"
	if got := runBed(t, "bed-resize", "-5", "100", "testdata/resize.bed"); got != want {
		t.Errorf("bed-resize mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestBedResizeRequiresExtent(t *testing.T) {
	// The bed-resize command and its flag vars are package globals reused
	// across Execute calls; reset the extent flags so this check does not
	// depend on test ordering.
	bedResizeLen5 = 0
	bedResizeLen3 = 0

	root := &cobra.Command{Use: "cgkit"}
	InitCmd(root)
	var buf bytes.Buffer
	root.SetOut(&buf)
	root.SetErr(&buf)
	root.SetArgs([]string{"bed-resize", "testdata/resize.bed"})
	if err := root.Execute(); err == nil {
		t.Fatal("expected an error when neither -5 nor -3 is set")
	}
}

func TestBedStats(t *testing.T) {
	want := "Total number of regions:\t7\n" +
		"Total number of bases:\t2260\n" +
		"\n" +
		"Mean region size:\t322.85714285714283\n" +
		"Median region size:\t100\n" +
		"Max region size:\t1000\n" +
		"Min region size:\t10\n" +
		"\n" +
		"chr1\t3\n" +
		"chr2\t1\n" +
		"chr3\t3\n"
	if got := runBed(t, "bed-stats", "testdata/stats.bed"); got != want {
		t.Errorf("bed-stats mismatch.\n got: %q\nwant: %q", got, want)
	}
}
