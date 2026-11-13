package vcfcmd

import (
	"os"
	"path/filepath"
	"reflect"
	"strings"
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
	if err := s.CatalogRuns(func(r varstore.CalledSiteRun) bool {
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

	vcfToVarstoreOut = out
	vcfToVarstoreMinDP = 10
	vcfToVarstoreBands = bands
	if len(bands) == 0 {
		vcfToVarstoreDepthRatio = 0 // unbanded
	}
	vcfToVarstoreForce = true
	t.Cleanup(func() {
		vcfToVarstoreOut, vcfToVarstoreMinDP, vcfToVarstoreBands, vcfToVarstoreForce = "", 10, nil, false
		vcfToVarstoreDepthRatio = varstore.DefaultDepthRatio
	})

	cmd := vcfToVarstoreCmd
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
	runs := runsOf(t, convertBands(t, []int{}))
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
	plain := runsOf(t, convertBands(t, []int{}))

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

// The worked example the reference-call gate was argued from. S1 is 0/0 at depth
// 25, 40 and 45 -- one catalog run in the 20-49 band with floor 25 -- and S2
// carries the alternate everywhere, so every position stays in the catalog.
const gateVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	1000	.	G	A	50	PASS	.	GT:DP	0/0:25	0/1:60
chr1	2000	.	G	A	50	PASS	.	GT:DP	0/0:40	0/1:60
chr1	3000	.	G	A	50	PASS	.	GT:DP	0/0:45	0/1:60
`

// homRefRows runs vcf-varquery --hom-ref over the whole contig at one gate and
// returns S1's reference rows as sample@pos. --sample rather than a
// region, because a region query on a plain VCF needs a tabix index.
func homRefRows(t *testing.T, input, gate string) map[string]bool {
	t.Helper()
	out := map[string]bool{}
	for _, l := range strings.Split(dataRowsOnly(runVcf(t, "vcf-varquery",
		"--sample", "S1", "--hom-ref", "--min-dp", gate, input)), "\n") {
		f := strings.Split(l, "\t")
		if len(f) == numCols && f[colChrom] != "chrom" && f[colGT] == varstore.HomRefGT {
			out[f[colSample]+"@"+f[colPos]] = true
		}
	}
	return out
}

// The backends' reference calls under a depth gate.
//
// A store keeps the run's floor, not each site's depth, so the promise is exact
// agreement wherever the gate cannot split a run -- at or below the conversion
// gate and on a band boundary -- and a subset everywhere else. A subset means the
// store may withhold a reference call the VCF makes; it must never make one the
// VCF does not, which is the direction that would fabricate an observation.
func TestHomRefBackendsAgreeOnBandBoundariesAndNestInside(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(in, []byte(gateVCF), 0o644); err != nil {
		t.Fatal(err)
	}
	store := convert(t, in, "--depth-bands", "10,20,50")

	for _, gate := range []string{"10", "20", "50"} {
		vcf, st := homRefRows(t, in, gate), homRefRows(t, store, gate)
		if !reflect.DeepEqual(vcf, st) {
			t.Errorf("gate %s is a boundary, so the backends must agree exactly:\n vcf   %v\n store %v",
				gate, vcf, st)
		}
	}

	// Inside the 20-49 band. The VCF knows S1 was 40x and 45x at two sites; the
	// store knows only that S1's run bottomed out at 25.
	vcf, st := homRefRows(t, in, "30"), homRefRows(t, store, "30")
	for row := range st {
		if !vcf[row] {
			t.Errorf("gate 30: the store reported %s as 0/0 and the VCF did not", row)
		}
	}
	if want := map[string]bool{"S1@2000": true, "S1@3000": true}; !reflect.DeepEqual(vcf, want) {
		t.Fatalf("fixture: the VCF's 0/0 rows at gate 30 = %v, want %v", vcf, want)
	}
	if len(st) != 0 {
		t.Errorf("gate 30: S1's run floor is 25, so the store should withhold every 0/0; got %v", st)
	}
}

// Generated bands are the default, and the two band flags are alternatives.
func TestDepthRatioIsTheDefaultAndExcludesExplicitBands(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(in, []byte(gateVCF), 0o644); err != nil {
		t.Fatal(err)
	}

	s, err := varstore.OpenParquet(convert(t, in))
	if err != nil {
		t.Fatal(err)
	}
	p := s.VolumeManifest().Params
	s.Close()
	want := []int32{10, 15, 22, 33, 50, 75, 113, 170, 256, 384, 576, 864}
	if !reflect.DeepEqual(p.DepthBands, want) || p.DepthRatio != varstore.DefaultDepthRatio {
		t.Errorf("default bands = %v ratio %v, want %v ratio %v",
			p.DepthBands, p.DepthRatio, want, varstore.DefaultDepthRatio)
	}

	s, err = varstore.OpenParquet(convert(t, in, "--depth-bands", "10,20,30"))
	if err != nil {
		t.Fatal(err)
	}
	p = s.VolumeManifest().Params
	s.Close()
	if !reflect.DeepEqual(p.DepthBands, []int32{10, 20, 30}) || p.DepthRatio != 0 {
		t.Errorf("explicit bands = %v ratio %v, want [10 20 30] ratio 0", p.DepthBands, p.DepthRatio)
	}

	err = runVcfErr(t, "vcf-tovarstore", "--out", filepath.Join(dir, "both"),
		"--depth-bands", "10,20", "--depth-ratio", "2", in)
	if err == nil || !strings.Contains(err.Error(), "alternatives") {
		t.Errorf("--depth-bands with --depth-ratio: err = %v, want a refusal", err)
	}
}
