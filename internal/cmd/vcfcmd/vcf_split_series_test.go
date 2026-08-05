package vcfcmd

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// A vcf-split series is written whole or not at all. Both guards below exist
// because vcf-concat --chunks stops at the first *missing* file rather than the
// first invalid one, so a truncated tail chunk or a leftover from a longer
// previous run is silently reassembled as though it were the entire input.

func TestSplitRefusesAnExistingSeries(t *testing.T) {
	base := filepath.Join(t.TempDir(), "out")
	runVcf(t, "vcf-split", "--out", base, "--num", "2", "testdata/sample.vcf")

	first, err := os.ReadFile(base + ".1.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}

	err = runVcfErr(t, "vcf-split", "--out", base, "--num", "2", "testdata/sample.vcf")
	if err == nil {
		t.Fatal("a second run overwrote an existing series without --force")
	}
	if !strings.Contains(err.Error(), "--force") {
		t.Errorf("error %q should say how to proceed", err)
	}
	// The refusal happens before anything is opened, so the previous series is
	// untouched rather than half-replaced.
	again, err := os.ReadFile(base + ".1.vcf.gz")
	if err != nil || string(again) != string(first) {
		t.Error("the existing series was modified by a refused run")
	}
}

// --force must clear the whole old series, not just the chunks the new run
// happens to overwrite. A shorter run would otherwise leave the old tail in
// place, and vcf-concat would read it as a continuation of the new output.
func TestForceRemovesStaleTailChunks(t *testing.T) {
	base := filepath.Join(t.TempDir(), "out")
	// --num 1 over a 5-record file gives 5 chunks.
	runVcf(t, "vcf-split", "--out", base, "--num", "1", "testdata/sample.vcf")
	if _, err := os.Stat(base + ".5.vcf.gz"); err != nil {
		t.Fatalf("expected 5 chunks from the first run: %v", err)
	}

	// --num 5 gives 1. Chunks 2..5 are now stale.
	runVcf(t, "vcf-split", "--force", "--out", base, "--num", "5", "testdata/sample.vcf")
	if _, err := os.Stat(base + ".1.vcf.gz"); err != nil {
		t.Fatalf("the new series is missing its first chunk: %v", err)
	}
	for n := 2; n <= 5; n++ {
		p := base + "." + string(rune('0'+n)) + ".vcf.gz"
		if _, err := os.Stat(p); err == nil {
			t.Errorf("stale chunk %s survived --force; vcf-concat would splice it in", p)
		}
	}
}

// A run that dies partway removes the chunks it wrote. A partial series is the
// dangerous outcome, not merely an untidy one.
func TestFailedSplitLeavesNoChunks(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "broken.vcf")
	const body = "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tG\t50\tPASS\t.\n" +
		"chr1\t200\t.\tC\tT\t50\tPASS\t.\n" +
		"THIS IS NOT A VCF RECORD\n"
	if err := os.WriteFile(in, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}

	base := filepath.Join(dir, "out")
	if err := runVcfErr(t, "vcf-split", "--out", base, "--num", "1", in); err == nil {
		t.Fatal("expected a parse error")
	}
	for n := 1; n <= 3; n++ {
		p := base + "." + string(rune('0'+n)) + ".vcf.gz"
		if _, err := os.Stat(p); err == nil {
			t.Errorf("a failed run left %s behind", p)
		}
	}
}

// A refused index does not abort the split. writeVcfTbi refuses to index an
// unsorted chunk after that chunk is already complete, and stopping there would
// leave exactly the partial series the guards above exist to prevent -- so the
// series is finished, indexes are skipped, and the error is reported at the end.
func TestUnsortedIndexStillFinishesTheSeries(t *testing.T) {
	base := filepath.Join(t.TempDir(), "out")
	err := runVcfErr(t, "vcf-split", "--out", base, "--num", "2", "--tbi", "testdata/unsorted.vcf")
	if err == nil {
		t.Fatal("expected an unsorted-index error")
	}
	// Every chunk the input calls for is present...
	for n := 1; n <= 2; n++ {
		p := base + "." + string(rune('0'+n)) + ".vcf.gz"
		if _, serr := os.Stat(p); serr != nil {
			t.Errorf("chunk %s is missing; the series was left partial: %v", p, serr)
		}
	}
	// ...and none of them was indexed.
	if _, serr := os.Stat(base + ".1.vcf.gz.tbi"); serr == nil {
		t.Error("an index was written for an unsorted chunk")
	}
}
