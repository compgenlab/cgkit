package vcfcmd

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/compgenlab/cghts/varstore"
)

// Banded callable runs.
//
// A run means "called at or above the gate throughout", and that alone makes a
// reference at 60x and one at exactly 10x indistinguishable forever -- which is
// fine until two sources disagree about a person and the reference side of the
// argument has nothing to say.
//
// Recording the lowest depth inside a run fixes that, and banding is what keeps
// the bound tight: without it, one poorly covered site drags a whole span down
// and every reference call across it inherits the worst moment in it.

// A sample climbing through the bands: 12, 15 (band 0), then 30, 35 (band 1),
// then 60 (band 2). Unbanded this is one run with MinDP 12; banded it is three,
// each bounding its own span.
const bandVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	G	A	50	PASS	.	GT:DP	0/0:12
chr1	200	.	G	A	50	PASS	.	GT:DP	0/0:15
chr1	300	.	G	A	50	PASS	.	GT:DP	0/0:30
chr1	400	.	G	A	50	PASS	.	GT:DP	0/0:35
chr1	500	.	G	A	50	PASS	.	GT:DP	0/0:60
`

func runsOf(t *testing.T, dir string) []varstore.CalledSiteRun {
	t.Helper()
	s, err := varstore.OpenParquet(dir)
	if err != nil {
		t.Fatalf("OpenParquet: %v", err)
	}
	defer s.Close()
	var out []varstore.CalledSiteRun
	if err := s.Regions(func(r varstore.CalledSiteRun) bool {
		out = append(out, r)
		return true
	}); err != nil {
		t.Fatalf("Regions: %v", err)
	}
	return out
}

func convertBands(t *testing.T, bands []int) string {
	t.Helper()
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(in, []byte(bandVCF), 0o644); err != nil {
		t.Fatal(err)
	}
	out := filepath.Join(dir, "store")

	vcfToParquetOut = out
	vcfToParquetMinDP = 10
	vcfToParquetBands = bands
	vcfToParquetForce = true
	t.Cleanup(func() {
		vcfToParquetOut, vcfToParquetMinDP, vcfToParquetBands, vcfToParquetForce = "", 10, []int{10, 20, 50}, false
	})

	cmd := vcfToParquetCmd
	cmd.SetArgs([]string{in})
	if err := cmd.RunE(cmd, []string{in}); err != nil {
		t.Fatalf("convert: %v", err)
	}
	return out
}

func TestBandedRunsBoundTheirOwnSpan(t *testing.T) {
	runs := runsOf(t, convertBands(t, []int{10, 20, 50}))
	if len(runs) != 3 {
		t.Fatalf("got %d runs, want 3 (one per depth band): %+v", len(runs), runs)
	}
	// Each run's MinDP must bound its own span rather than the whole sample.
	want := []struct{ start, end, minDP int32 }{
		{100, 200, 12},
		{300, 400, 30},
		{500, 500, 60},
	}
	for i, w := range want {
		r := runs[i]
		if r.Start != w.start || r.End != w.end || r.MinDP != w.minDP {
			t.Errorf("run %d = start %d end %d minDP %d, want %d %d %d",
				i, r.Start, r.End, r.MinDP, w.start, w.end, w.minDP)
		}
	}
}

// Unbanded, the same input is one run whose bound is the worst site in it --
// which is exactly the looseness banding exists to avoid, and is still correct.
func TestUnbandedRunsTakeTheWorstDepth(t *testing.T) {
	runs := runsOf(t, convertBands(t, nil))
	if len(runs) != 1 {
		t.Fatalf("got %d runs, want 1 when unbanded: %+v", len(runs), runs)
	}
	if runs[0].MinDP != 12 {
		t.Errorf("MinDP = %d, want 12 -- the lowest depth anywhere in the run", runs[0].MinDP)
	}
	if runs[0].Start != 100 || runs[0].End != 500 {
		t.Errorf("run spans %d-%d, want 100-500", runs[0].Start, runs[0].End)
	}
}

// Banding must not change WHO is callable, only how the coverage is described.
// A boundary that quietly dropped sites would move every count on the store.
func TestBandingPreservesCallableSites(t *testing.T) {
	banded := runsOf(t, convertBands(t, []int{10, 20, 50}))
	plain := runsOf(t, convertBands(t, nil))

	var bandedSites, plainSites int32
	for _, r := range banded {
		bandedSites += r.NSites
	}
	for _, r := range plain {
		plainSites += r.NSites
	}
	if bandedSites != plainSites {
		t.Errorf("banded covers %d sites, unbanded %d -- banding must only describe "+
			"coverage differently, never change it", bandedSites, plainSites)
	}
}
