package vcfcmd

import (
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

// homRefSamplesAt runs vcf-varquery --hom-ref for one exact locus and returns the
// samples reported 0/0.
func homRefSamplesAt(t *testing.T, input, locus string) []string {
	t.Helper()
	var out []string
	for _, l := range strings.Split(dataRowsOnly(runVcf(t, "vcf-varquery", "--hom-ref",
		"--min-dp", "10", "--variant", locus, input)), "\n") {
		f := strings.Split(l, "\t")
		if len(f) == numCols && f[colChrom] != "chrom" && f[colGT] == "0/0" {
			out = append(out, f[colSample])
		}
	}
	return out
}

// THE END-TO-END CHECK the series was built toward. A store converted from the
// pVCF, with each sample's gVCF extracted and attached, must answer an
// off-catalog position the way the gVCFs themselves do when queried directly.
// Fragments are --passthrough here, so the store holds exactly the source's
// blocks and the two must agree sample for sample; coarsening may legitimately
// bridge a short gap a gVCF leaves open, which is a different, recorded claim.
func TestAttachedStoreAnswersOffCatalogLikeTheGvcfs(t *testing.T) {
	roster, s1 := refBlocksFixture(t, "S1")
	dir := t.TempDir()
	// S2 is S1's layout with the 300-399 block at 3x: below the gate.
	s2 := filepath.Join(dir, "S2.g.vcf")
	g2 := strings.Replace(refBlocksGVCF("S2"),
		"END=399\tGT:DP:GQ:MIN_DP\t0/0:30:30:30", "END=399\tGT:DP:GQ:MIN_DP\t0/0:3:30:3", 1)
	if err := os.WriteFile(s2, []byte(g2), 0o644); err != nil {
		t.Fatal(err)
	}

	frags := filepath.Join(dir, "frags")
	for name, g := range map[string]string{"S1": s1, "S2": s2} {
		runVcf(t, "vcf-refblocks", "--roster", roster, "--passthrough",
			"--out", filepath.Join(frags, name), g)
	}
	runVcf(t, "vcf-varstore", "attach-refblocks", roster, frags)

	gvcfs := map[string]string{}
	for name, g := range map[string]string{"S1": s1, "S2": s2} {
		gz := filepath.Join(dir, name+".g.vcf.gz")
		runVcf(t, "vcf-filter", "-o", gz, "--tbi", g)
		gvcfs[name] = gz
	}

	for _, locus := range []string{
		"chr1:150:C:T",  // a block both share
		"chr1:350:C:T",  // S2's block there is below the gate
		"chr1:602:C:T",  // inside a 3x block: nobody
		"chr1:650:C:T",  // a block both share
		"chr1:705:C:T",  // just past S1's own SNV at 700
		"chr1:5000:C:T", // no block at all
	} {
		var direct []string
		for _, name := range []string{"S1", "S2"} {
			direct = append(direct, homRefSamplesAt(t, gvcfs[name], locus)...)
		}
		store := homRefSamplesAt(t, roster, locus)
		if !reflect.DeepEqual(store, direct) {
			t.Errorf("%s: store reports %v as 0/0, the gVCFs %v", locus, store, direct)
		}
	}
}

// A chromosome with no catalog site belongs to no volume, so blocks on it have
// nowhere to go: the fixture's catalog is chr1 only, and S1's chr2 blocks are not
// attached. Taking them anyway would put every chromosome's blocks into every
// volume of an archive, whose pVCF headers all declare every contig.
func TestBlocksOnAChromosomeNoVolumeHoldsAreNotAttached(t *testing.T) {
	roster, s1 := refBlocksFixture(t, "S1")
	frag := filepath.Join(t.TempDir(), "S1")
	runVcf(t, "vcf-refblocks", "--roster", roster, "--out", frag, s1)
	runVcf(t, "vcf-varstore", "attach-refblocks", roster, frag)
	if got := homRefSamplesAt(t, roster, "chr2:150:C:T"); len(got) != 0 {
		t.Errorf("chr2 is in no volume, yet the store reports %v as 0/0 there", got)
	}
}

func TestAttachRefusesASecondAttachWithoutRebuildOrReplace(t *testing.T) {
	roster, s1 := refBlocksFixture(t, "S1")
	frag := filepath.Join(t.TempDir(), "S1")
	runVcf(t, "vcf-refblocks", "--roster", roster, "--out", frag, s1)
	runVcf(t, "vcf-varstore", "attach-refblocks", roster, frag)

	if err := runVcfErr(t, "vcf-varstore", "attach-refblocks", roster, frag); err == nil ||
		!strings.Contains(err.Error(), "already has ref_blocks") {
		t.Errorf("a second attach: err = %v", err)
	}
	runVcf(t, "vcf-varstore", "attach-refblocks", "--rebuild", roster, frag)
	runVcf(t, "vcf-varstore", "attach-refblocks", "--replace", "S1", roster, frag)

	// Removing S1 leaves an empty table, so nobody is reference off the catalog.
	runVcf(t, "vcf-varstore", "attach-refblocks", "--replace", "S1", roster)
	if got := homRefSamplesAt(t, roster, "chr1:150:C:T"); len(got) != 0 {
		t.Errorf("after removing S1's blocks, chr1:150 still reports %v", got)
	}
}

func TestAttachExpandsADirectoryAndAList(t *testing.T) {
	dir := t.TempDir()
	for _, n := range []string{"a", "b"} {
		if err := os.MkdirAll(filepath.Join(dir, "frags", n), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(filepath.Join(dir, "frags", n, "fragment.json"), []byte("{}"), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	list := filepath.Join(dir, "list.txt")
	if err := os.WriteFile(list, []byte("# fragments\nx/one\n\nx/two\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	got, err := expandFragments([]string{filepath.Join(dir, "frags")}, list)
	if err != nil {
		t.Fatal(err)
	}
	want := []string{filepath.Join(dir, "frags", "a"), filepath.Join(dir, "frags", "b"), "x/one", "x/two"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %v, want %v", got, want)
	}
	if _, err := expandFragments([]string{t.TempDir()}, ""); err == nil {
		t.Error("an empty directory was accepted as a source of fragments")
	}
}

func TestVarsummaryReportsAttachedRefBlocks(t *testing.T) {
	roster, s1 := refBlocksFixture(t, "S1")
	frag := filepath.Join(t.TempDir(), "S1")
	runVcf(t, "vcf-refblocks", "--roster", roster, "--out", frag, s1)
	runVcf(t, "vcf-varstore", "attach-refblocks", roster, frag)
	out := runVcf(t, "vcf-varsummary", roster)
	for _, want := range []string{"ref blocks yes", "gaps up to 10 bp", "ref_blocks"} {
		if !strings.Contains(out, want) {
			t.Errorf("vcf-varsummary is missing %q:\n%s", want, out)
		}
	}
}
