package samcmd

import (
	"io"
	"testing"

	"github.com/compgenlab/cghts/htsio"
)

func loadStats(t *testing.T, path string, opt statsOpts) *statsResult {
	t.Helper()
	reader, err := htsio.NewSamReader(path)
	if err != nil {
		t.Fatalf("NewSamReader(%s): %v", path, err)
	}
	defer reader.Close()
	header, err := reader.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	res, err := computeStats(reader, header, opt, io.Discard)
	if err != nil {
		t.Fatalf("computeStats: %v", err)
	}
	return res
}

func TestSamStatsBasicCounts(t *testing.T) {
	res := loadStats(t, "testdata/test.bam", statsOpts{})

	if res.total != 5 {
		t.Errorf("total = %d, want 5", res.total)
	}
	if res.mapped != 5 {
		t.Errorf("mapped = %d, want 5", res.mapped)
	}
	if res.unmapped != 0 {
		t.Errorf("unmapped = %d, want 0", res.unmapped)
	}
	if res.multiple != 0 {
		t.Errorf("multiple = %d, want 0", res.multiple)
	}
	if res.totalBases != 250 {
		t.Errorf("totalBases = %d, want 250 (5 reads x 50M)", res.totalBases)
	}
	if res.refLength != 150000 {
		t.Errorf("refLength = %d, want 150000 (chr1 100000 + chr2 50000)", res.refLength)
	}
	if res.refCounts["chr1"] != 4 || res.refCounts["chr2"] != 1 {
		t.Errorf("refCounts = chr1:%d chr2:%d, want chr1:4 chr2:1", res.refCounts["chr1"], res.refCounts["chr2"])
	}
}

func TestSamStatsSecondOfPairAndUnmapped(t *testing.T) {
	// test_tags.bam has 11 records; read_paired_2 (flag 0x80) is skipped as the
	// second read of a pair, and read_unmapped (flag 0x4) is counted but unmapped.
	res := loadStats(t, "testdata/test_tags.bam", statsOpts{})

	if res.total != 10 {
		t.Errorf("total = %d, want 10 (11 records minus second-of-pair)", res.total)
	}
	if res.mapped != 9 {
		t.Errorf("mapped = %d, want 9", res.mapped)
	}
	if res.unmapped != 1 {
		t.Errorf("unmapped = %d, want 1", res.unmapped)
	}
	// 10 mapped 10M reads, all with phred-40 ('I') bases → 100% Q30.
	if res.totalBases != 100 {
		t.Errorf("totalBases = %d, want 100", res.totalBases)
	}
	if res.q30Bases != 100 {
		t.Errorf("q30Bases = %d, want 100 (all bases are phred 40)", res.q30Bases)
	}
}

func TestSamStatsTagTally(t *testing.T) {
	opt := statsOpts{tags: []tagSpec{{name: "NM", numeric: true}}}
	res := loadStats(t, "testdata/test_tags.bam", opt)

	nm := res.tagCounts["NM"]
	for _, val := range []string{"0", "1", "2", "3"} {
		if nm[val] != 1 {
			t.Errorf("NM[%s] = %d, want 1", val, nm[val])
		}
	}
	// 5 mapped reads reaching the tag step lack an NM tag.
	if res.tagMissing["NM"] != 5 {
		t.Errorf("NM missing = %d, want 5", res.tagMissing["NM"])
	}
}

func TestSamStatsRgidFilter(t *testing.T) {
	res := loadStats(t, "testdata/test_tags.bam", statsOpts{rgid: "sample1"})

	// sample1 reads: multi_tags, paired_1, paired_2 (skipped, 2nd of pair),
	// chr2, unmapped → total 4, of which one is unmapped.
	if res.total != 4 {
		t.Errorf("total = %d, want 4 for rgid=sample1", res.total)
	}
	if res.mapped != 3 {
		t.Errorf("mapped = %d, want 3 for rgid=sample1", res.mapped)
	}
	if res.unmapped != 1 {
		t.Errorf("unmapped = %d, want 1 for rgid=sample1", res.unmapped)
	}
}

func TestParseTagSpecs(t *testing.T) {
	specs, err := parseTagSpecs("NH:i,RG:Z,MAPQ")
	if err != nil {
		t.Fatalf("parseTagSpecs: %v", err)
	}
	want := []tagSpec{{"NH", true}, {"RG", false}, {"MAPQ", true}}
	if len(specs) != len(want) {
		t.Fatalf("len(specs) = %d, want %d", len(specs), len(want))
	}
	for i := range want {
		if specs[i] != want[i] {
			t.Errorf("specs[%d] = %+v, want %+v", i, specs[i], want[i])
		}
	}

	if _, err := parseTagSpecs("NH"); err == nil {
		t.Error("expected error for tag without type")
	}
	if _, err := parseTagSpecs("XF:f"); err == nil {
		t.Error("expected error for unsupported tag type")
	}
}
