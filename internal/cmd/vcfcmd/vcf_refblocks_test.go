package vcfcmd

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/varstore"
)

// A pVCF with S1 and S2, converted with bands 10,20,50: the store a fragment is
// built for.
const refBlocksPVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##contig=<ID=chr2,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	700	.	A	G	50	PASS	.	GT:DP	0/1:40	0/0:30
`

// S1's gVCF, built to hit each rule once. Expected spans, at --min-dp 10, bands
// 10,20,50 and --max-gap 10:
//
//	chr1 100-299 @12   two band-10 blocks that abut
//	chr1 300-499 @30   a band change opens it; a gap of exactly 10 is bridged
//	chr1 511-699 @40   a gap of 11 is not; a 3x block inside is dropped and bridged
//	chr1 701-799 @40   the sample's own SNV at 700 is never bridged
//	chr2 100-199 @25   a second chromosome
func refBlocksGVCF(sample string) string {
	block := func(chrom string, start, end, dp int) string {
		return fmt.Sprintf("%s\t%d\t.\tA\t<NON_REF>\t.\t.\tEND=%d\tGT:DP:GQ:MIN_DP\t0/0:%d:30:%d\n",
			chrom, start, end, dp, dp)
	}
	return `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##contig=<ID=chr2,length=1000000>
##ALT=<ID=NON_REF,Description="Any other allele">
##INFO=<ID=END,Number=1,Type=Integer,Description="End of the block">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth in the block">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	` + sample + "\n" +
		block("chr1", 100, 199, 12) +
		block("chr1", 200, 299, 15) +
		block("chr1", 300, 399, 30) +
		block("chr1", 410, 499, 35) +
		block("chr1", 511, 599, 40) +
		block("chr1", 600, 604, 3) +
		block("chr1", 605, 699, 40) +
		"chr1\t700\t.\tA\tG,<NON_REF>\t50\tPASS\t.\tGT:DP:GQ\t0/1:40:99\n" +
		block("chr1", 701, 799, 40) +
		block("chr2", 100, 199, 25)
}

func refBlocksFixture(t *testing.T, sample string) (roster, gvcf string) {
	t.Helper()
	dir := t.TempDir()
	pvcf := filepath.Join(dir, "cohort.vcf")
	gvcf = filepath.Join(dir, sample+".g.vcf")
	if err := os.WriteFile(pvcf, []byte(refBlocksPVCF), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(gvcf, []byte(refBlocksGVCF(sample)), 0o644); err != nil {
		t.Fatal(err)
	}
	return convert(t, pvcf, "--depth-bands", "10,20,50"), gvcf
}

func fragmentSpans(t *testing.T, dir string) (*varstore.FragmentMeta, []string) {
	t.Helper()
	m, err := varstore.ReadFragmentMeta(context.Background(), dir)
	if err != nil {
		t.Fatal(err)
	}
	var spans []string
	for _, c := range m.Chroms {
		if err := varstore.FragmentBlocks(context.Background(), dir, c, func(b varstore.RefBlock) bool {
			spans = append(spans, fmt.Sprintf("%s %d-%d @%d", b.Chrom, b.Start, b.End, b.MinDP))
			return true
		}); err != nil {
			t.Fatal(err)
		}
	}
	return m, spans
}

func TestRefBlocksCoarsensAGvcf(t *testing.T) {
	roster, gvcf := refBlocksFixture(t, "S1")
	out := filepath.Join(t.TempDir(), "S1.frag")
	runVcf(t, "vcf-refblocks", "--roster", roster, "--out", out, gvcf)

	m, spans := fragmentSpans(t, out)
	want := []string{
		"chr1 100-299 @12",
		"chr1 300-499 @30",
		"chr1 511-699 @40",
		"chr1 701-799 @40",
		"chr2 100-199 @25",
	}
	if !reflect.DeepEqual(spans, want) {
		t.Errorf("spans:\n got %v\nwant %v", spans, want)
	}
	// Settings come from the store unless named, so the fragment can be attached.
	if m.SampleID != "S1" || m.MinDP != 10 || !reflect.DeepEqual(m.DepthBands, []int32{10, 20, 50}) || m.MaxGap != 10 {
		t.Errorf("meta = %+v, want S1 at the store's min-dp 10 and bands 10,20,50, max-gap 10", m)
	}
}

func TestRefBlocksPassthroughKeepsSourceBlocks(t *testing.T) {
	roster, gvcf := refBlocksFixture(t, "S1")
	out := filepath.Join(t.TempDir(), "S1.frag")
	runVcf(t, "vcf-refblocks", "--roster", roster, "--out", out, "--passthrough", gvcf)
	_, spans := fragmentSpans(t, out)
	// Every source block above the gate, unmerged: 9 blocks less the 3x one.
	if len(spans) != 8 {
		t.Errorf("passthrough kept %d spans, want 8: %v", len(spans), spans)
	}
}

func TestRefBlocksRefusals(t *testing.T) {
	roster, gvcf := refBlocksFixture(t, "S1")

	_, offRoster := refBlocksFixture(t, "S9")
	if err := runVcfErr(t, "vcf-refblocks", "--roster", roster,
		"--out", filepath.Join(t.TempDir(), "x"), offRoster); err == nil ||
		!strings.Contains(err.Error(), "not among") {
		t.Errorf("a sample off the roster: err = %v", err)
	}

	if err := runVcfErr(t, "vcf-refblocks", "--roster", roster, "--min-dp", "20",
		"--out", filepath.Join(t.TempDir(), "x"), gvcf); err == nil ||
		!strings.Contains(err.Error(), "disagrees") {
		t.Errorf("a --min-dp the store was not built at: err = %v", err)
	}

	if err := runVcfErr(t, "vcf-refblocks", "--roster", roster, "--depth-bands", "10,30",
		"--out", filepath.Join(t.TempDir(), "x"), gvcf); err == nil ||
		!strings.Contains(err.Error(), "disagree") {
		t.Errorf("bands the store was not built with: err = %v", err)
	}

	if err := runVcfErr(t, "vcf-refblocks", "--roster", gvcf,
		"--out", filepath.Join(t.TempDir(), "x"), gvcf); err == nil {
		t.Error("a VCF accepted as the roster store")
	}

	out := filepath.Join(t.TempDir(), "S1.frag")
	runVcf(t, "vcf-refblocks", "--roster", roster, "--out", out, gvcf)
	if err := runVcfErr(t, "vcf-refblocks", "--roster", roster, "--out", out, gvcf); err == nil ||
		!strings.Contains(err.Error(), "--force") {
		t.Errorf("rewriting a finished fragment without --force: err = %v", err)
	}
	runVcf(t, "vcf-refblocks", "--roster", roster, "--out", out, "--force", gvcf)
}
