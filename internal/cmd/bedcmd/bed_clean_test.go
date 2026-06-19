package bedcmd

import (
	"bytes"
	"testing"

	"github.com/spf13/cobra"
)

// golden output captured from `ngsutilsj bed-clean testdata/input.bed`.
const bedCleanGolden = "chr1\t100\t200\t\t0\t+\n" +
	"chr2\t50\t75\tregionB\t2\t-\n" +
	"chr3\t10\t20\tregionC\t0\t+\n" +
	"chr1\t100\t200\tregionA\t1\t+\tfoo\tbar\n"

func TestBedClean(t *testing.T) {
	root := &cobra.Command{Use: "cgio"}
	InitCmd(root)

	var buf bytes.Buffer
	root.SetOut(&buf)
	root.SetErr(&buf)
	root.SetArgs([]string{"bed-clean", "testdata/input.bed"})

	if err := root.Execute(); err != nil {
		t.Fatalf("Execute: %v", err)
	}

	if got := buf.String(); got != bedCleanGolden {
		t.Errorf("bed-clean output mismatch.\n got: %q\nwant: %q", got, bedCleanGolden)
	}
}
