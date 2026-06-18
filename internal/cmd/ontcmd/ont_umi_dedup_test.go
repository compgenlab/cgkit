package ontcmd

import (
	"os"
	"strings"
	"testing"

	"github.com/compgenlab/hts/htsio"
	"github.com/compgenlab/hts/htsio/bam"
)

// makeTestBAM writes a coordinate-sorted BAM to path with the given records.
func makeTestBAM(t *testing.T, path string, records []*htsio.SamRecord) {
	t.Helper()
	header := htsio.NewSamHeader()
	header.AddLine("@HD\tVN:1.6\tSO:coordinate")
	header.AddLine("@SQ\tSN:chr1\tLN:100000")
	header.AddLine("@SQ\tSN:chr2\tLN:100000")

	writer, err := bam.NewWriter(path, header)
	if err != nil {
		t.Fatalf("NewSamWriter: %v", err)
	}
	for _, rec := range records {
		if err := writer.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close writer: %v", err)
	}
}

// readAllBAM reads all records from a BAM file.
func readAllBAM(t *testing.T, path string) []*htsio.SamRecord {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	reader, err := bam.NewReader(f, path, nil)
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

	var recs []*htsio.SamRecord
	for rec, err := range reader.Records() {
		if err != nil {
			t.Fatalf("Records: %v", err)
		}
		recs = append(recs, rec)
	}
	return recs
}

// rec is a helper to build a SamRecord with common fields. SEQ is sized to the
// CIGAR's query length so the record is well-formed (the BAM/CRAM writers reject
// a CIGAR/SEQ length mismatch); seq supplies the base pattern, repeated or
// truncated to the required length.
func rec(name string, flag int, ref string, pos int, cigar string, seq string, tags map[string]htsio.SamTag) *htsio.SamRecord {
	if tags == nil {
		tags = make(map[string]htsio.SamTag)
	}
	return &htsio.SamRecord{
		ReadName:  name,
		Flag:      flag,
		RefName:   ref,
		Pos:       pos,
		MapQ:      60,
		Cigar:     cigar,
		RefNext:   "*",
		PosNext:   0,
		InsertLen: 0,
		Seq:       fitSeq(seq, htsio.CigarQueryLen(cigar)),
		Qual:      "*",
		Tags:      tags,
	}
}

// fitSeq returns a sequence of exactly n bases built by repeating pattern
// (truncating the final copy). It returns "*" for n == 0 (no query bases).
func fitSeq(pattern string, n int) string {
	if n == 0 {
		return "*"
	}
	if pattern == "" || pattern == "*" {
		pattern = "A"
	}
	var b strings.Builder
	for b.Len() < n {
		b.WriteString(pattern)
	}
	return b.String()[:n]
}

func tags(kv ...interface{}) map[string]htsio.SamTag {
	m := make(map[string]htsio.SamTag)
	for i := 0; i < len(kv); i += 3 {
		key := kv[i].(string)
		typ := kv[i+1].(int32) // Go character literals are int32 (rune)
		val := kv[i+2].(string)
		m[key] = htsio.SamTag{Type: byte(typ), Value: val}
	}
	return m
}

// --- Unit tests for selector logic ---

func TestParseTagSelectorFlag(t *testing.T) {
	tests := []struct {
		input     string
		wantTag   string
		wantAsc   bool
		wantError bool
	}{
		{"AS+", "AS", false, false},
		{"AS", "AS", false, false},  // no suffix defaults to higher wins
		{"NM-", "NM", true, false},
		{"XY+", "XY", false, false},
		{"A", "", false, true},      // too short
	}
	for _, tt := range tests {
		ts, err := parseTagSelectorFlag(tt.input)
		if tt.wantError {
			if err == nil {
				t.Errorf("parseTagSelectorFlag(%q) expected error", tt.input)
			}
			continue
		}
		if err != nil {
			t.Errorf("parseTagSelectorFlag(%q) unexpected error: %v", tt.input, err)
			continue
		}
		if ts.tag != tt.wantTag || ts.ascending != tt.wantAsc {
			t.Errorf("parseTagSelectorFlag(%q) = {tag:%q, asc:%v}, want {tag:%q, asc:%v}",
				tt.input, ts.tag, ts.ascending, tt.wantTag, tt.wantAsc)
		}
	}
}

func TestSelectBest_SingleTagHigherWins(t *testing.T) {
	// AS+ means higher alignment score wins.
	sel := []selector{&tagSelector{tag: "AS", ascending: false}}

	reads := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "50M", "ACGT", tags("AS", 'i', "100", "MI", 'Z', "mi_1")),
		rec("r2", 0, "chr1", 100, "50M", "ACGT", tags("AS", 'i', "200", "MI", 'Z', "mi_1")),
		rec("r3", 0, "chr1", 100, "50M", "ACGT", tags("AS", 'i', "150", "MI", 'Z', "mi_1")),
	}
	best := selectBest(reads, sel)
	if reads[best].ReadName != "r2" {
		t.Errorf("selectBest AS+ = %q, want r2 (AS=200)", reads[best].ReadName)
	}
}

func TestSelectBest_SingleTagLowerWins(t *testing.T) {
	// NM- means lower edit distance wins.
	sel := []selector{&tagSelector{tag: "NM", ascending: true}}

	reads := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "50M", "ACGT", tags("NM", 'i', "5", "MI", 'Z', "mi_1")),
		rec("r2", 0, "chr1", 100, "50M", "ACGT", tags("NM", 'i', "1", "MI", 'Z', "mi_1")),
		rec("r3", 0, "chr1", 100, "50M", "ACGT", tags("NM", 'i', "3", "MI", 'Z', "mi_1")),
	}
	best := selectBest(reads, sel)
	if reads[best].ReadName != "r2" {
		t.Errorf("selectBest NM- = %q, want r2 (NM=1)", reads[best].ReadName)
	}
}

func TestSelectBest_ChainedSelectors(t *testing.T) {
	// AS+ first, then NM- to break ties.
	sel := []selector{
		&tagSelector{tag: "AS", ascending: false},
		&tagSelector{tag: "NM", ascending: true},
	}

	reads := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "50M", "ACGT", tags("AS", 'i', "200", "NM", 'i', "5", "MI", 'Z', "mi_1")),
		rec("r2", 0, "chr1", 100, "50M", "ACGT", tags("AS", 'i', "200", "NM", 'i', "2", "MI", 'Z', "mi_1")),
		rec("r3", 0, "chr1", 100, "50M", "ACGT", tags("AS", 'i', "100", "NM", 'i', "0", "MI", 'Z', "mi_1")),
	}
	best := selectBest(reads, sel)
	// r1 and r2 tied on AS=200, r2 wins on NM=2 < NM=5.
	if reads[best].ReadName != "r2" {
		t.Errorf("selectBest AS+,NM- = %q, want r2 (AS=200,NM=2)", reads[best].ReadName)
	}
}

func TestSelectBest_LongestTiebreaker(t *testing.T) {
	// AS+ first, then longest (aligned query length, excluding soft clips).
	sel := []selector{
		&tagSelector{tag: "AS", ascending: false},
		&longestSelector{},
	}

	reads := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "50M", "ACGTACGT", tags("AS", 'i', "200", "MI", 'Z', "mi_1")),
		rec("r2", 0, "chr1", 100, "80M", "ACGTACGTACGT", tags("AS", 'i', "200", "MI", 'Z', "mi_1")),
	}
	best := selectBest(reads, sel)
	// Tied on AS=200, r2 wins because longer aligned query (80M > 50M).
	if reads[best].ReadName != "r2" {
		t.Errorf("selectBest AS+,longest = %q, want r2 (longer aligned)", reads[best].ReadName)
	}
}

func TestSelectBest_LongestExcludesSoftClip(t *testing.T) {
	sel := []selector{&longestSelector{}}

	reads := []*htsio.SamRecord{
		// r1: 20S + 80M = 80 aligned query bases (longer SEQ but less aligned)
		rec("r1", 0, "chr1", 100, "20S80M", "ACGTACGTACGTACGTACGT", tags("MI", 'Z', "mi_1")),
		// r2: 90M = 90 aligned query bases
		rec("r2", 0, "chr1", 100, "90M", "ACGTACGTACGT", tags("MI", 'Z', "mi_1")),
	}
	best := selectBest(reads, sel)
	// r2 has more aligned bases (90 vs 80) despite shorter SEQ.
	if reads[best].ReadName != "r2" {
		t.Errorf("selectBest longest = %q, want r2 (90M > 20S80M)", reads[best].ReadName)
	}
}

func TestSelectBest_MissingTag(t *testing.T) {
	// Read with missing tag should lose to reads that have it.
	sel := []selector{&tagSelector{tag: "AS", ascending: false}}

	reads := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "50M", "ACGT", tags("MI", 'Z', "mi_1")),                       // no AS
		rec("r2", 0, "chr1", 100, "50M", "ACGT", tags("AS", 'i', "50", "MI", 'Z', "mi_1")),      // AS=50
	}
	best := selectBest(reads, sel)
	if reads[best].ReadName != "r2" {
		t.Errorf("selectBest with missing tag = %q, want r2", reads[best].ReadName)
	}
}

// --- Integration tests using BAM round-trip ---

func TestUmiDedup_BasicSelection(t *testing.T) {
	dir := t.TempDir()
	inPath := dir + "/input.bam"
	outPath := dir + "/output.bam"

	// Three reads in MI group "mi_000001.001", one read in "mi_000001.002".
	// In group 1: r2 has the highest AS (200).
	// Group 2: only r4, so it should always be kept.
	input := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_000001.001", "AS", 'i', "100", "NM", 'i', "5")),
		rec("r2", 0, "chr1", 110, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_000001.001", "AS", 'i', "200", "NM", 'i', "2")),
		rec("r3", 0, "chr1", 120, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_000001.001", "AS", 'i', "150", "NM", 'i', "3")),
		rec("r4", 0, "chr1", 5000, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_000001.002", "AS", 'i', "80", "NM", 'i', "10")),
	}
	makeTestBAM(t, inPath, input)

	// Select by AS+ (highest wins).
	umiDedupOutput = outPath
	umiDedupMITag = "MI"
	umiDedupMarkDuplicates = false
	selectors := []selector{&tagSelector{tag: "AS", ascending: false}}

	if err := runUmiDedup(inPath, selectors, nil); err != nil {
		t.Fatalf("runUmiDedup: %v", err)
	}

	recs := readAllBAM(t, outPath)
	names := make(map[string]bool)
	for _, r := range recs {
		names[r.ReadName] = true
	}

	if !names["r2"] {
		t.Error("expected r2 (AS=200) to be kept")
	}
	if !names["r4"] {
		t.Error("expected r4 (sole member of group 2) to be kept")
	}
	if names["r1"] || names["r3"] {
		t.Errorf("expected r1, r3 to be removed; got names: %v", names)
	}
}

func TestUmiDedup_MarkDuplicates(t *testing.T) {
	dir := t.TempDir()
	inPath := dir + "/input.bam"
	outPath := dir + "/output.bam"

	input := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "100")),
		rec("r2", 0, "chr1", 110, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "200")),
	}
	makeTestBAM(t, inPath, input)

	umiDedupOutput = outPath
	umiDedupMITag = "MI"
	umiDedupMarkDuplicates = true
	selectors := []selector{&tagSelector{tag: "AS", ascending: false}}

	if err := runUmiDedup(inPath, selectors, nil); err != nil {
		t.Fatalf("runUmiDedup: %v", err)
	}

	recs := readAllBAM(t, outPath)
	if len(recs) != 2 {
		t.Fatalf("expected 2 records with --mark-duplicates, got %d", len(recs))
	}

	for _, r := range recs {
		isDup := r.Flag&0x400 != 0
		if r.ReadName == "r2" && isDup {
			t.Error("r2 (best) should NOT have dup flag")
		}
		if r.ReadName == "r1" && !isDup {
			t.Error("r1 (non-best) should have dup flag")
		}
	}
}

func TestUmiDedup_ChainedSelectorsIntegration(t *testing.T) {
	dir := t.TempDir()
	inPath := dir + "/input.bam"
	outPath := dir + "/output.bam"

	// r1 and r2 tied on AS=200. r2 has lower NM (1 < 3), so r2 wins.
	input := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "200", "NM", 'i', "3")),
		rec("r2", 0, "chr1", 110, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "200", "NM", 'i', "1")),
	}
	makeTestBAM(t, inPath, input)

	umiDedupOutput = outPath
	umiDedupMITag = "MI"
	umiDedupMarkDuplicates = false
	selectors := []selector{
		&tagSelector{tag: "AS", ascending: false},
		&tagSelector{tag: "NM", ascending: true},
	}

	if err := runUmiDedup(inPath, selectors, nil); err != nil {
		t.Fatalf("runUmiDedup: %v", err)
	}

	recs := readAllBAM(t, outPath)
	if len(recs) != 1 {
		t.Fatalf("expected 1 record, got %d", len(recs))
	}
	if recs[0].ReadName != "r2" {
		t.Errorf("expected r2, got %q", recs[0].ReadName)
	}
}

func TestUmiDedup_NoMIPassthrough(t *testing.T) {
	dir := t.TempDir()
	inPath := dir + "/input.bam"
	outPath := dir + "/output.bam"

	// r_noMI has no MI tag and should pass through. r1 is the sole MI group member.
	input := []*htsio.SamRecord{
		rec("r_noMI", 0, "chr1", 50, "50M", "ACGT", tags("AS", 'i', "50")),
		rec("r1", 0, "chr1", 100, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "100")),
	}
	makeTestBAM(t, inPath, input)

	umiDedupOutput = outPath
	umiDedupMITag = "MI"
	umiDedupMarkDuplicates = false
	selectors := []selector{&tagSelector{tag: "AS", ascending: false}}

	if err := runUmiDedup(inPath, selectors, nil); err != nil {
		t.Fatalf("runUmiDedup: %v", err)
	}

	recs := readAllBAM(t, outPath)
	names := make(map[string]bool)
	for _, r := range recs {
		names[r.ReadName] = true
	}
	if !names["r_noMI"] {
		t.Error("read without MI should pass through")
	}
	if !names["r1"] {
		t.Error("sole MI group member should be kept")
	}
}

func TestUmiDedup_SecondarySupplementaryDropped(t *testing.T) {
	dir := t.TempDir()
	inPath := dir + "/input.bam"
	outPath := dir + "/output.bam"

	// All secondary/supplementary reads are dropped. Only primaries are kept.
	input := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "100")),
		rec("r2", 0, "chr1", 110, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "200")),
		rec("r1", 0x100, "chr1", 150, "50M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "50")),    // secondary for r1
		rec("r2", 0x800, "chr1", 160, "60M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "80")),    // supplementary for r2
	}
	makeTestBAM(t, inPath, input)

	umiDedupOutput = outPath
	umiDedupMITag = "MI"
	umiDedupMarkDuplicates = false
	selectors := []selector{&tagSelector{tag: "AS", ascending: false}}

	if err := runUmiDedup(inPath, selectors, nil); err != nil {
		t.Fatalf("runUmiDedup: %v", err)
	}

	recs := readAllBAM(t, outPath)
	// Should have only r2 primary — all sec/supp dropped.
	if len(recs) != 1 {
		var names []string
		for _, r := range recs {
			names = append(names, r.ReadName)
		}
		t.Fatalf("expected 1 record (r2 primary only), got %d: %v", len(recs), names)
	}
	if recs[0].ReadName != "r2" {
		t.Errorf("expected r2, got %q", recs[0].ReadName)
	}
}

func TestUmiDedup_MultipleChromosomes(t *testing.T) {
	dir := t.TempDir()
	inPath := dir + "/input.bam"
	outPath := dir + "/output.bam"

	// Two MI groups on different chromosomes, same MI value (shouldn't happen
	// in practice from ont-umi-cluster, but tests chromosome flush logic).
	input := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "100")),
		rec("r2", 0, "chr1", 110, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "200")),
		rec("r3", 0, "chr2", 100, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_2", "AS", 'i', "300")),
		rec("r4", 0, "chr2", 110, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_2", "AS", 'i', "150")),
	}
	makeTestBAM(t, inPath, input)

	umiDedupOutput = outPath
	umiDedupMITag = "MI"
	umiDedupMarkDuplicates = false
	selectors := []selector{&tagSelector{tag: "AS", ascending: false}}

	if err := runUmiDedup(inPath, selectors, nil); err != nil {
		t.Fatalf("runUmiDedup: %v", err)
	}

	recs := readAllBAM(t, outPath)
	names := make(map[string]bool)
	for _, r := range recs {
		names[r.ReadName] = true
	}
	if !names["r2"] {
		t.Error("expected r2 (best in chr1 MI group)")
	}
	if !names["r3"] {
		t.Error("expected r3 (best in chr2 MI group)")
	}
	if len(recs) != 2 {
		t.Errorf("expected 2 records, got %d", len(recs))
	}
}

func TestUmiDedup_LongestOnly(t *testing.T) {
	dir := t.TempDir()
	inPath := dir + "/input.bam"
	outPath := dir + "/output.bam"

	input := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "50M", "ACGT", tags("MI", 'Z', "mi_1")),
		rec("r2", 0, "chr1", 110, "120M", "ACGTACGTACGT", tags("MI", 'Z', "mi_1")),
		rec("r3", 0, "chr1", 120, "80M", "ACGTACGT", tags("MI", 'Z', "mi_1")),
	}
	makeTestBAM(t, inPath, input)

	umiDedupOutput = outPath
	umiDedupMITag = "MI"
	umiDedupMarkDuplicates = false
	selectors := []selector{&longestSelector{}}

	if err := runUmiDedup(inPath, selectors, nil); err != nil {
		t.Fatalf("runUmiDedup: %v", err)
	}

	recs := readAllBAM(t, outPath)
	if len(recs) != 1 {
		t.Fatalf("expected 1 record, got %d", len(recs))
	}
	if recs[0].ReadName != "r2" {
		t.Errorf("expected r2 (longest), got %q", recs[0].ReadName)
	}
}

func TestUmiDedup_StatsReport(t *testing.T) {
	dir := t.TempDir()
	inPath := dir + "/input.bam"
	outPath := dir + "/output.bam"
	statsPath := dir + "/stats.txt"

	// 3 reads in group mi_1 (AS: 100, 200, 150; NM: 5, 2, 3).
	// 1 read in group mi_2.
	// 1 read without MI (passthrough).
	input := []*htsio.SamRecord{
		rec("r_noMI", 0, "chr1", 50, "50M", "ACGT", tags("AS", 'i', "50")),
		rec("r1", 0, "chr1", 100, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "100", "NM", 'i', "5")),
		rec("r2", 0, "chr1", 110, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "200", "NM", 'i', "2")),
		rec("r3", 0, "chr1", 120, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_1", "AS", 'i', "150", "NM", 'i', "3")),
		rec("r4", 0, "chr1", 5000, "100M", "ACGTACGTAC", tags("MI", 'Z', "mi_2", "AS", 'i', "80", "NM", 'i', "10")),
	}
	makeTestBAM(t, inPath, input)

	umiDedupOutput = outPath
	umiDedupMITag = "MI"
	umiDedupMarkDuplicates = false
	umiDedupStatsFile = statsPath
	defer func() { umiDedupStatsFile = "" }()

	selectors := []selector{
		&tagSelector{tag: "AS", ascending: false},
		&tagSelector{tag: "NM", ascending: true},
	}

	if err := runUmiDedup(inPath, selectors, []string{"AS", "NM"}); err != nil {
		t.Fatalf("runUmiDedup: %v", err)
	}

	// Verify stats file was written and contains expected content.
	data, err := os.ReadFile(statsPath)
	if err != nil {
		t.Fatalf("reading stats file: %v", err)
	}
	content := string(data)

	// Check key lines are present.
	checks := []string{
		"Total reads:        5",
		"Primary:          4",
		"No MI (passthru): 1",
		"MI groups:          2",
		"Reads kept:         2",
		"Reads discarded:    2",
		"Duplication rate:   50.0%",
		"AS tag distribution",
		"NM tag distribution",
		"kept",
		"discarded",
	}
	for _, check := range checks {
		if !contains(content, check) {
			t.Errorf("stats file missing %q\nfull content:\n%s", check, content)
		}
	}

	// Check group size histogram: should have one group of size 3, one of size 1.
	if !contains(content, "1\t1") || !contains(content, "3\t1") {
		t.Errorf("stats file missing expected group size histogram entries\nfull content:\n%s", content)
	}
}

func contains(s, substr string) bool {
	return len(s) >= len(substr) && searchString(s, substr)
}

func searchString(s, substr string) bool {
	for i := 0; i <= len(s)-len(substr); i++ {
		if s[i:i+len(substr)] == substr {
			return true
		}
	}
	return false
}
