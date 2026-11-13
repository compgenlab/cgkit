package vcfcmd

import (
	"context"
	"github.com/compgenlab/cghts/varstore"
	"os"
	"path/filepath"
	"testing"
)

// writeFile drops a fixture into a temp dir and returns its path.
func writeFile(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

// TestSiteListWorksForBothCommands pins that one file serves both commands.
//
// They now share a single parser rather than two that had to agree, so this checks
// the remaining seam: vcf-gtcount's conversion of the shared target set into its
// own site type.
func TestSiteListWorksForBothCommands(t *testing.T) {
	list := writeFile(t, "shared.txt",
		"# shared site list\nchr1 100 A G\nchr1 500 A T\nchr2:300:G:C\n")

	targets, err := varstore.ParseTargets(context.Background(), []string{list})
	if err != nil {
		t.Fatalf("target parser: %v", err)
	}
	if len(targets.Loci) != 3 {
		t.Fatalf("want 3 loci from the shared list, got %d", len(targets.Loci))
	}

	sites, err := collectGtSites(context.Background(), nil, []string{list})
	if err != nil {
		t.Fatalf("vcf-gtcount rejected the shared list: %v", err)
	}
	if len(sites) != len(targets.Loci) {
		t.Fatalf("ParseTargets read %d targets, vcf-gtcount %d", len(targets.Loci), len(sites))
	}
	for i, l := range targets.Loci {
		g := sites[i]
		if l.Chrom != g.chrom || int(l.Pos) != g.pos || l.Ref != g.ref || l.Alt != g.alt {
			t.Errorf("row %d differs: ParseTargets %s, gtcount %s:%d:%s:%s",
				i, l, g.chrom, g.pos, g.ref, g.alt)
		}
	}
}
