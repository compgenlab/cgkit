package tabcmd

import (
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/compgenlab/cghts/htsio/bgzf"
	"github.com/compgenlab/cghts/htsio/tabix"
	"github.com/spf13/cobra"
)

func runTab(t *testing.T, args ...string) {
	t.Helper()
	root := &cobra.Command{Use: "cgkit"}
	InitCmd(root)
	root.SetArgs(args)
	if err := root.Execute(); err != nil {
		t.Fatalf("Execute(%v): %v", args, err)
	}
}

func TestTabixIndexCommand(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "regions.bed.gz")

	// Write a sorted, BGZF-compressed BED with no index.
	w, err := bgzf.NewBGZipFile(path)
	if err != nil {
		t.Fatal(err)
	}
	for _, l := range []string{"chr1\t90\t110\tgeneA", "chr2\t400\t600\tgeneC"} {
		io.WriteString(w, l+"\n")
	}
	w.Close()

	// Index it via the command.
	runTab(t, "tabix-index", "--preset", "bed", path)

	if _, err := os.Stat(path + ".tbi"); err != nil {
		t.Fatalf("expected %s.tbi to be created: %v", path, err)
	}

	// The index must support random access (including the single chr2 record).
	r, err := tabix.NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	seq, err := r.Query("chr2", 499, 500)
	if err != nil {
		t.Fatal(err)
	}
	n := 0
	for rec, err := range seq {
		if err != nil {
			t.Fatal(err)
		}
		if rec.Line != "chr2\t400\t600\tgeneC" {
			t.Errorf("got %q", rec.Line)
		}
		n++
	}
	if n != 1 {
		t.Errorf("chr2 query returned %d rows, want 1", n)
	}
}
