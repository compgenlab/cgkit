package vcfcmd

import (
	"strings"
	"testing"
)

func TestVcfToBedpe(t *testing.T) {
	out := runVcf(t, "vcf-tobedpe", "testdata/sv.vcf")
	if !strings.HasPrefix(out, "#chrom1\tstart1\tstop1\tchrom2\tstart2\tstop2\tname\n") {
		t.Errorf("missing BEDPE header:\n%s", out)
	}
	// Deletion: start1=POS, stop1=POS+1, partner offset back by one.
	if !strings.Contains(out, "chr1\t50\t51\tchr1\t119\t120\t<DEL>\n") {
		t.Errorf("DEL row wrong:\n%s", out)
	}
	// BND: point interval at POS-1..POS, partner from the bracket.
	if !strings.Contains(out, "chr1\t99\t100\tchr1\t199\t200\tN[chr1:200[\n") {
		t.Errorf("BND row wrong:\n%s", out)
	}
}

func TestVcfToBedpeNameScore(t *testing.T) {
	out := runVcf(t, "vcf-tobedpe", "--name", "@ID", "--score", "SVTYPE:INFO", "testdata/sv.vcf")
	if !strings.Contains(out, "\tname\tscore\n") {
		t.Errorf("--score should add a score column:\n%s", out)
	}
	// name = ID (sv_del), score = INFO SVTYPE (DEL).
	if !strings.Contains(out, "\tsv_del\tDEL\n") {
		t.Errorf("--name @ID / --score wrong:\n%s", out)
	}
}

func TestVcfSvToFasta(t *testing.T) {
	out := runVcf(t, "vcf-svtofasta", "--bnd", "--flanking", "10", "testdata/svref.fa", "testdata/sv.vcf")
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	// Four BND records -> four >sv records (the <DEL> is skipped).
	var headers, seqs int
	for _, l := range lines {
		if strings.HasPrefix(l, ">sv|") {
			headers++
		} else if l != "" {
			seqs++
		}
	}
	if headers != 4 || seqs != 4 {
		t.Errorf("expected 4 BND records, got %d headers / %d seqs:\n%s", headers, seqs, out)
	}
	if !strings.Contains(out, ">sv|chr1|100|chr1|200|bnd_a|BND|3to5\n") {
		t.Errorf("missing bnd_a header:\n%s", out)
	}
	// flanking=10 -> two 10bp flanks joined into a 20bp sequence.
	if !strings.Contains(out, "\nCTAAACCTATAAACCTTTCT\n") {
		t.Errorf("bnd_a sequence wrong:\n%s", out)
	}
}

func TestVcfSvToFastaNeedsIndex(t *testing.T) {
	// Missing .fai should be a clear error (sample.vcf has no companion .fa).
	if err := runVcfErr(t, "vcf-svtofasta", "--bnd", "testdata/missing.fa", "testdata/sv.vcf"); err == nil {
		t.Errorf("expected an error for a missing FAI index")
	}
}
