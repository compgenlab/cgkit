package fastqcmd

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/compgenlab/cgkit/internal/cmdtest"
)

func TestFastqGC(t *testing.T) {
	p := filepath.Join(t.TempDir(), "in.fq")
	const fq = "@r1\nGGCC\n+\nIIII\n" + // all G/C
		"@r2\nAATT\n+\nIIII\n" + // none
		"@r3\nACGT\n+\nIIII\n" // half
	if err := os.WriteFile(p, []byte(fq), 0o644); err != nil {
		t.Fatal(err)
	}
	got := cmdtest.Run(t, InitCmd, "fastq-gc", p)
	want := "r1\t1.0000\nr2\t0.0000\nr3\t0.5000\n"
	if got != want {
		t.Errorf("fastq-gc mismatch.\n got: %q\nwant: %q", got, want)
	}
}

// fastqcmd.InitCmd used to rely on fastacmd.InitCmd having registered the
// shared help group first, so calling it alone panicked inside cobra -- which
// is why the one test this package had built its root by hand.
func TestInitCmdStandsAlone(t *testing.T) {
	root, _ := cmdtest.NewRoot(InitCmd)
	if !root.ContainsGroup("fastaqcmd") {
		t.Fatal("fastqcmd.InitCmd did not register its own help group")
	}
	// The real check: cobra panics when a command names a group the root does
	// not have, and it does so while rendering help.
	root.SetArgs([]string{"--help"})
	if err := root.Execute(); err != nil {
		t.Fatalf("help failed for a fastqcmd-only root: %v", err)
	}
}
