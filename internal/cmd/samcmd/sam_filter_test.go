package samcmd

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/hts/htsio"
	"github.com/spf13/cobra"
)

// TestSamFilterProvenance verifies that sam-filter stamps a @PG provenance line
// (with PN, VN, and CL) into the output header. The output is written as BAM
// (the default SAM-text output goes to stdout), then the header is read back.
func TestSamFilterProvenance(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.sam")
	out := filepath.Join(dir, "out.bam")
	const sam = "@HD\tVN:1.6\n" +
		"@SQ\tSN:chr1\tLN:1000\n" +
		"r1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\t*\n"
	if err := os.WriteFile(in, []byte(sam), 0o644); err != nil {
		t.Fatal(err)
	}

	root := &cobra.Command{Use: "cgkit"}
	InitCmd(root)
	root.SetArgs([]string{"sam-filter", "--bam", in, out})
	if err := root.Execute(); err != nil {
		t.Fatalf("sam-filter: %v", err)
	}

	reader, err := htsio.NewSamReader(out)
	if err != nil {
		t.Fatalf("NewSamReader(%s): %v", out, err)
	}
	defer reader.Close()
	header, err := reader.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}

	var pg string
	for _, line := range header.Lines {
		if strings.HasPrefix(line, "@PG\t") && strings.Contains(line, "ID:sam-filter") {
			pg = line
		}
	}
	if pg == "" {
		t.Fatalf("no sam-filter @PG provenance line in output header:\n%v", header.Lines)
	}
	for _, want := range []string{"PN:cgkit", "VN:", "CL:"} {
		if !strings.Contains(pg, want) {
			t.Errorf("@PG line missing %q: %s", want, pg)
		}
	}
}
