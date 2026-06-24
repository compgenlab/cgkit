package vcfcmd

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestVcfFilterValue exercises the batch-2 value/flag filters: numeric INFO
// comparisons, an INFO flag-present check, and a per-sample string equality.
func TestVcfFilterValue(t *testing.T) {
	// DP is an INFO field (30,10,25,25,40). --lt DP:25:INFO flags only DP<25.
	// The "INFO" sample name is appended to the generated ID (as in ngsutilsj),
	// and a flagged row is rewritten so QUAL loses its trailing ".0".
	out := runVcf(t, "vcf-filter", "--lt", "DP:25:INFO", "testdata/annotate.vcf")
	if !strings.Contains(out, "chr1\t140\t.\tC\tT\t20\tlowqual;DP_lt_25.0_INFO\t") {
		t.Errorf("chr1:140 should be flagged DP_lt_25.0_INFO:\n%s", out)
	}
	if !strings.Contains(out, "chr1\t100\t.\tA\tG\t50.0\tPASS\t") {
		t.Errorf("chr1:100 (DP=30) should stay PASS:\n%s", out)
	}
	if !strings.Contains(out, "##FILTER=<ID=DP_lt_25.0_INFO,") {
		t.Errorf("missing DP_lt_25.0_INFO FILTER def:\n%s", out)
	}

	// --gte DP:30:INFO flags DP>=30 (chr1:100 and chr2:500).
	out = runVcf(t, "vcf-filter", "--gte", "DP:30:INFO", "testdata/annotate.vcf")
	if !strings.Contains(out, "chr1\t100\t.\tA\tG\t50\tDP_gte_30.0_INFO\t") {
		t.Errorf("chr1:100 should be flagged DP_gte_30.0_INFO:\n%s", out)
	}
	if !strings.Contains(out, "chr2\t500\t.\tCAT\tC\t99\tDP_gte_30.0_INFO\t") {
		t.Errorf("chr2:500 should be flagged DP_gte_30.0_INFO:\n%s", out)
	}
	if !strings.Contains(out, "chr1\t150\t.\tA\tC\t40.0\tPASS\t") {
		t.Errorf("chr1:150 (DP=25) should stay PASS:\n%s", out)
	}

	// DP is present in every record's INFO, so --flag-present DP flags them all.
	out = runVcf(t, "vcf-filter", "--flag-present", "DP", "testdata/annotate.vcf")
	if strings.Contains(out, "\tPASS\t") {
		t.Errorf("--flag-present DP should flag every record:\n%s", out)
	}
	if !strings.Contains(out, "##FILTER=<ID=DP_present,") {
		t.Errorf("missing DP_present FILTER def:\n%s", out)
	}

	// TUMOR GT is 0/1 except chr1:150 (1/1). --eq GT:0/1:TUMOR flags the 0/1 rows.
	out = runVcf(t, "vcf-filter", "--eq", "GT:0/1:TUMOR", "testdata/annotate.vcf")
	if !strings.Contains(out, "chr1\t100\t.\tA\tG\t50\tGT_eq_0/1_TUMOR\t") {
		t.Errorf("chr1:100 (TUMOR 0/1) should be flagged:\n%s", out)
	}
	if !strings.Contains(out, "chr1\t150\t.\tA\tC\t40.0\tPASS\t") {
		t.Errorf("chr1:150 (TUMOR 1/1) should stay PASS:\n%s", out)
	}
}

// TestVcfFilterStats checks the --stats per-combination tally file.
func TestVcfFilterStats(t *testing.T) {
	statsFile := filepath.Join(t.TempDir(), "stats.txt")
	runVcf(t, "vcf-filter", "--snv", "--qual", "30", "--chrom-fail", "chr2",
		"--stats", statsFile, "testdata/annotate.vcf")
	data, err := os.ReadFile(statsFile)
	if err != nil {
		t.Fatalf("reading stats file: %v", err)
	}
	got := string(data)
	want := "SNV\t2\n" +
		"QUAL_lt_30.0,SNV,lowqual\t1\n" +
		"CHROM_FAIL_chr2\t1\n"
	if got != want {
		t.Errorf("stats mismatch.\n got: %q\nwant: %q", got, want)
	}
}
