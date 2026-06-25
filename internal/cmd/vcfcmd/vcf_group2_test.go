package vcfcmd

import (
	"path/filepath"
	"strings"
	"testing"
)

func TestVcfConcat(t *testing.T) {
	out := runVcf(t, "vcf-concat", "testdata/concat_a.vcf", "testdata/concat_b.vcf")
	var pos []string
	for _, line := range strings.Split(strings.TrimRight(out, "\n"), "\n") {
		if strings.HasPrefix(line, "#") {
			continue
		}
		f := strings.Split(line, "\t")
		pos = append(pos, f[0]+":"+f[1])
	}
	want := "chr1:100,chr1:200,chr1:300,chr1:400,chr2:500,chr2:600"
	if strings.Join(pos, ",") != want {
		t.Errorf("concat order = %v, want %s", pos, want)
	}
}

func TestVcfConcatOverlap(t *testing.T) {
	// Concatenating a file with itself duplicates every position -> error.
	if err := runVcfErr(t, "vcf-concat", "testdata/concat_a.vcf", "testdata/concat_a.vcf"); err == nil {
		t.Errorf("expected an overlapping-position error")
	} else if !strings.Contains(err.Error(), "overlapping") {
		t.Errorf("error = %v, want overlapping-position", err)
	}
}

func TestVcfConcatChunks(t *testing.T) {
	// Split sample.vcf into numbered chunks, then recombine with --chunks.
	base := filepath.Join(t.TempDir(), "split")
	runVcf(t, "vcf-split", "--out", base, "--num", "2", "testdata/sample.vcf")

	out := runVcf(t, "vcf-concat", "--chunks", base+".1.vcf.gz")
	var pos []string
	for _, line := range strings.Split(strings.TrimRight(out, "\n"), "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		f := strings.Split(line, "\t")
		pos = append(pos, f[0]+":"+f[1])
	}
	want := "chr1:100,chr1:200,chr1:300,chr2:500,chr2:1000"
	if strings.Join(pos, ",") != want {
		t.Errorf("--chunks order = %v, want %s", pos, want)
	}
}

func TestVcfMerge(t *testing.T) {
	out := runVcf(t, "vcf-merge", "testdata/merge_a.vcf", "testdata/merge_b.vcf")
	// Primary ID kept; INFO and FORMAT unioned (primary first).
	if !strings.Contains(out, "chr1\t100\trs100\tA\tG\t50\tPASS\tAA=1;BB=2\tGT:AD:DP\t0/0:28,2:30\t0/1:15,15:30\n") {
		t.Errorf("merged chr1:100 wrong:\n%s", out)
	}
	// The union'd header carries both annotation defs.
	if !strings.Contains(out, "##INFO=<ID=AA,") || !strings.Contains(out, "##INFO=<ID=BB,") ||
		!strings.Contains(out, "##FORMAT=<ID=AD,") || !strings.Contains(out, "##FORMAT=<ID=DP,") {
		t.Errorf("merged header missing unioned defs:\n%s", out)
	}
}

// TestVcfMergeTakesID covers cgio's ID-from-secondary merge (the intended
// behavior; ngsutilsj NPEs whenever the primary ID is missing).
func TestVcfMergeTakesID(t *testing.T) {
	out := runVcf(t, "vcf-merge", "testdata/merge_noid.vcf", "testdata/merge_b.vcf")
	if !strings.Contains(out, "chr1\t100\trs1\tA\tG\t") || !strings.Contains(out, "chr1\t200\trs2\tC\tT\t") {
		t.Errorf("missing primary ID should be taken from the secondary:\n%s", out)
	}
}

func TestVcfMergeOutOfOrder(t *testing.T) {
	// concat_b has different variants than merge_a -> lockstep mismatch.
	if err := runVcfErr(t, "vcf-merge", "testdata/merge_a.vcf", "testdata/concat_b.vcf"); err == nil {
		t.Errorf("expected a variants-out-of-order error")
	} else if !strings.Contains(err.Error(), "out of order") && !strings.Contains(err.Error(), "extra sample") {
		t.Errorf("error = %v", err)
	}
}
