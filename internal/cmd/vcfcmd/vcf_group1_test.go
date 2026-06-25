package vcfcmd

import (
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/compgenlab/hts/vcf"
)

func TestVcfClearFilter(t *testing.T) {
	out := runVcf(t, "vcf-clearfilter", "--filter", "lowqual", "testdata/sample.vcf")
	if !strings.Contains(out, "chr1\t200\t.\tC\tT\t20\tPASS\t") {
		t.Errorf("chr1:200 should clear to PASS:\n%s", out)
	}
	if !strings.Contains(out, "CG_CLEARED_FILTER=lowqual") {
		t.Errorf("CG_CLEARED_FILTER not recorded:\n%s", out)
	}
	if !strings.Contains(out, "##INFO=<ID=CG_CLEARED_FILTER,") {
		t.Errorf("missing CG_CLEARED_FILTER def:\n%s", out)
	}
}

func TestVcfRename(t *testing.T) {
	out := runVcf(t, "vcf-rename", "--sample", "NORMAL:GERMLINE", "--sample", "2:CANCER", "testdata/sample.vcf")
	if !strings.Contains(out, "\tGERMLINE\tCANCER\n") {
		t.Errorf("sample columns not renamed:\n%s", firstLines(out, 14))
	}
}

func TestVcfChrFix(t *testing.T) {
	out := runVcf(t, "vcf-chrfix", "--ensembl", "testdata/sample.vcf")
	if !strings.Contains(out, "\n1\t100\t") || strings.Contains(out, "\nchr1\t") {
		t.Errorf("--ensembl should drop the chr prefix:\n%s", out)
	}
	if !strings.Contains(out, "##contig=<ID=1,") {
		t.Errorf("contig header not converted:\n%s", out)
	}
	// --contigs keeps only the listed (post-conversion) chroms.
	out = runVcf(t, "vcf-chrfix", "--ensembl", "--contigs", "1", "testdata/sample.vcf")
	if strings.Contains(out, "\n2\t") {
		t.Errorf("--contigs 1 should drop chr2:\n%s", out)
	}
}

func TestVcfRemoveFlags(t *testing.T) {
	out := runVcf(t, "vcf-remove-flags", "testdata/sample.vcf")
	if !strings.Contains(out, "FLAGS=DB") {
		t.Errorf("DB flag not converted to FLAGS:\n%s", out)
	}
	if strings.Contains(out, ";DB\t") || strings.Contains(out, "=DB;") {
		t.Errorf("original DB flag still present:\n%s", out)
	}
	if !strings.Contains(out, "##INFO=<ID=FLAGS,") || strings.Contains(out, "##INFO=<ID=DB,") {
		t.Errorf("header defs not rewritten:\n%s", out)
	}
}

func TestVcfHeaderInfo(t *testing.T) {
	if got := runVcf(t, "vcf-header-info", "--filters", "testdata/sample.vcf"); got != "PASS\tAll filters passed\nlowqual\tLow quality, score < 30\n" {
		t.Errorf("--filters mismatch:\n%q", got)
	}
	if got := runVcf(t, "vcf-header-info", "--sample", "testdata/sample.vcf"); got != "NORMAL\nTUMOR\n" {
		t.Errorf("--sample mismatch:\n%q", got)
	}
	if got := runVcf(t, "vcf-header-info", "--contig", "testdata/sample.vcf"); got != "chr1\t248956422\nchr2\t242193529\n" {
		t.Errorf("--contig mismatch:\n%q", got)
	}
}

func TestVcfCheck(t *testing.T) {
	// runVcf merges stdout+stderr; the only output is the stderr summary.
	if got := runVcf(t, "vcf-check", "testdata/sample.vcf"); got != "OK: 5 variants\n" {
		t.Errorf("vcf-check output = %q, want the OK summary", got)
	}
}

func TestVcfSplit(t *testing.T) {
	base := filepath.Join(t.TempDir(), "out")
	runVcf(t, "vcf-split", "--out", base, "--num", "2", "testdata/sample.vcf")
	var recs int
	for i := 1; ; i++ {
		name := base + "." + strconv.Itoa(i) + ".vcf.gz"
		if _, err := os.Stat(name); err != nil {
			if i == 1 {
				t.Fatalf("no output files produced")
			}
			if i-1 != 3 {
				t.Errorf("expected 3 chunks for 5 vars/2, got %d", i-1)
			}
			break
		}
		recs += countVcfRecords(t, name)
	}
	if recs != 5 {
		t.Errorf("total records across chunks = %d, want 5", recs)
	}
}

func TestVcfSampleExport(t *testing.T) {
	out := runVcf(t, "vcf-sample-export", "--gt", "--id", "testdata/sample.vcf")
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if lines[0] != "chrom\tpos\tID\tref\talt\tsample\tGT" {
		t.Errorf("header = %q", lines[0])
	}
	if lines[1] != "chr1\t100\trs1\tA\tG\tNORMAL\tA/A" {
		t.Errorf("row 1 = %q", lines[1])
	}
	if lines[2] != "chr1\t100\trs1\tA\tG\tTUMOR\tA/G" {
		t.Errorf("row 2 = %q", lines[2])
	}
}

func TestVcfToCount(t *testing.T) {
	out := runVcf(t, "vcf-tocount", "--af", "--total", "testdata/annotate.vcf")
	if !strings.Contains(out, "chrom\tpos\tref\talt\tref_count\talt_count\talt_freq\ttotal_count\n") {
		t.Errorf("header missing:\n%s", out)
	}
	// chr1:100 NORMAL AD=14,1 -> ref 14, alt 1, af 1/15, total 15
	if !strings.Contains(out, "chr1\t100\tA\tG\t14\t1\t"+javaDouble(1.0/15.0)+"\t15\n") {
		t.Errorf("count row wrong:\n%s", out)
	}
}

func TestVcfStripAll(t *testing.T) {
	out := runVcf(t, "vcf-strip", "--all", "testdata/sample.vcf")
	if !strings.Contains(out, "chr1\t100\t.\tA\tG\t50\tPASS\t.\n") {
		t.Errorf("--all should leave 8 bare columns:\n%s", out)
	}
	if strings.Contains(out, "GT") || strings.Contains(out, "##INFO") || strings.Contains(out, "##FORMAT") {
		t.Errorf("--all left annotations behind:\n%s", out)
	}
	if strings.Contains(dataChrom(out), "FORMAT") {
		t.Errorf("#CHROM line still has sample columns:\n%s", out)
	}
}

func TestVcfStripKeep(t *testing.T) {
	// Remove all INFO except DP.
	out := runVcf(t, "vcf-strip", "--info", "*", "--keep-info", "DP", "testdata/sample.vcf")
	if !strings.Contains(out, "DP=30") {
		t.Errorf("DP should be kept:\n%s", out)
	}
	if strings.Contains(out, "AF=") || strings.Contains(out, ";DB") {
		t.Errorf("other INFO should be stripped:\n%s", out)
	}
	if !strings.Contains(out, "##INFO=<ID=DP,") || strings.Contains(out, "##INFO=<ID=AF,") {
		t.Errorf("INFO header defs not stripped correctly:\n%s", out)
	}
}

// --- helpers ---

func firstLines(s string, n int) string {
	lines := strings.Split(s, "\n")
	if len(lines) > n {
		lines = lines[:n]
	}
	return strings.Join(lines, "\n")
}

func dataChrom(out string) string {
	for _, line := range strings.Split(out, "\n") {
		if strings.HasPrefix(line, "#CHROM") {
			return line
		}
	}
	return ""
}

func countVcfRecords(t *testing.T, filename string) int {
	t.Helper()
	r, err := vcf.NewVcfFile(filename)
	if err != nil {
		t.Fatalf("open %s: %v", filename, err)
	}
	defer r.Close()
	if _, err := r.Header(); err != nil {
		t.Fatalf("header %s: %v", filename, err)
	}
	n := 0
	for {
		_, err := r.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("read %s: %v", filename, err)
		}
		n++
	}
	return n
}
