package ontcmd

// Tests for ont-polya. Three things to know about the shared BAM test harness in
// ont_umi_dedup_test.go (makeTestBAM, rec, fitSeq, tags), which these reuse:
//   - rec() runs seq through fitSeq, which *repeats* the pattern to fill the
//     CIGAR's query length, so pass exact-length strings to get a specific SEQ.
//   - The BAM writer rejects len(Seq) != CigarQueryLen(Cigar).
//   - The polya* flag vars are global mutable state; polyaTestDefaults restores
//     them via t.Cleanup.

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio"
)

// polyaNonA returns n bases containing neither A nor T, so it works as filler on
// either strand without perturbing the tail search.
func polyaNonA(n int) string {
	if n <= 0 {
		return ""
	}
	return strings.Repeat("CG", (n+1)/2)[:n]
}

func polyaTestParams() polyaParams {
	return polyaParams{window: 10, minAFrac: 0.8, minLen: 10, maxJunk: 20}
}

func TestFindPolyA(t *testing.T) {
	p := polyaTestParams()

	cases := []struct {
		name      string
		seq       string
		mRNAPlus  bool
		wantFound bool
		wantStart int
		wantEnd   int
		wantLen   int
	}{
		{
			name:      "sense clean tail",
			seq:       polyaNonA(30) + strings.Repeat("A", 20),
			mRNAPlus:  true,
			wantFound: true, wantStart: 30, wantEnd: 49, wantLen: 20,
		},
		{
			name:      "sense tail with one error",
			seq:       polyaNonA(30) + "AAAAAGAAAAAAAAAA",
			mRNAPlus:  true,
			wantFound: true, wantStart: 30, wantEnd: 45, wantLen: 16,
		},
		{
			// A scattered N is bridged by the window, same as any other miscall.
			name:      "sense tail with N",
			seq:       polyaNonA(30) + "AAAAANAAAAAAAAAA",
			mRNAPlus:  true,
			wantFound: true, wantStart: 30, wantEnd: 45, wantLen: 16,
		},
		{
			// An N run abutting the tail must NOT extend it: N counts against the
			// A-fraction, so the trace stops at index 25. This is the case that
			// catches sequtils.DNAMatches creeping in — it treats N as matching A
			// (ConvertDNATo4Bit('N') is 0x0F), which would drag the start back to
			// 20 and report length 15. A scattered N cannot catch that, since the
			// window bridges it either way.
			name:      "sense N run before tail",
			seq:       polyaNonA(20) + strings.Repeat("N", 5) + strings.Repeat("A", 10),
			mRNAPlus:  true,
			wantFound: true, wantStart: 25, wantEnd: 34, wantLen: 10,
		},
		{
			// A miscall near the 3' tip must not cost the anchor. The window there
			// covers only a few real bases, so one T four bases from the end is
			// 3/4 = 0.75 and fails -- unless the denominator stays pinned at
			// `window`. With a partial denominator the anchor is lost, the retry
			// re-anchors inward, and this reports length 16 instead of 20.
			name:      "sense noisy tail tip",
			seq:       polyaNonA(50) + "AAAGAAAAACAAAAAATAAA",
			mRNAPlus:  true,
			wantFound: true, wantStart: 50, wantEnd: 69, wantLen: 20,
		},
		{
			name:     "sense no tail",
			seq:      strings.Repeat("ACGT", 12),
			mRNAPlus: true,
		},
		{
			name:     "sense short tail",
			seq:      polyaNonA(40) + strings.Repeat("A", 5),
			mRNAPlus: true,
		},
		{
			name:      "sense exactly minLen",
			seq:       polyaNonA(40) + strings.Repeat("A", 10),
			mRNAPlus:  true,
			wantFound: true, wantStart: 40, wantEnd: 49, wantLen: 10,
		},
		{
			// The retry case: without it the scan anchors on the adapter's lone
			// trailing A, fails the window on the C, and reports no call.
			name:      "sense junk then tail",
			seq:       polyaNonA(30) + strings.Repeat("A", 20) + "GGATCA",
			mRNAPlus:  true,
			wantFound: true, wantStart: 30, wantEnd: 49, wantLen: 20,
		},
		{
			name:     "sense junk beyond budget",
			seq:      strings.Repeat("A", 20) + polyaNonA(30),
			mRNAPlus: true,
		},
		{
			name:      "sense all A",
			seq:       strings.Repeat("A", 50),
			mRNAPlus:  true,
			wantFound: true, wantStart: 0, wantEnd: 49, wantLen: 50,
		},
		{
			name:      "sense lowercase tail",
			seq:       polyaNonA(30) + strings.Repeat("a", 10),
			mRNAPlus:  true,
			wantFound: true, wantStart: 30, wantEnd: 39, wantLen: 10,
		},
		{
			// The window at index 20 is exactly 8/10 A and index 20 is itself an
			// A, so the tail extends to 20 only if the threshold is inclusive.
			// A ">" comparison would stop at 23 and report length 8.
			name:      "sense window exactly at threshold",
			seq:       polyaNonA(20) + "A" + "CG" + strings.Repeat("A", 8),
			mRNAPlus:  true,
			wantFound: true, wantStart: 20, wantEnd: 30, wantLen: 11,
		},
		{
			name:      "antisense clean tail",
			seq:       strings.Repeat("T", 20) + polyaNonA(30),
			mRNAPlus:  false,
			wantFound: true, wantStart: 19, wantEnd: 0, wantLen: 20,
		},
		{
			name:      "antisense junk then tail",
			seq:       "TGATCC" + strings.Repeat("T", 20) + polyaNonA(30),
			mRNAPlus:  false,
			wantFound: true, wantStart: 25, wantEnd: 6, wantLen: 20,
		},
		{
			name:      "antisense all T",
			seq:       strings.Repeat("T", 50),
			mRNAPlus:  false,
			wantFound: true, wantStart: 49, wantEnd: 0, wantLen: 50,
		},
		{
			name:     "empty seq",
			seq:      "",
			mRNAPlus: true,
		},
		{
			name:     "unspecified seq",
			seq:      "*",
			mRNAPlus: true,
		},
		{
			name:     "seq shorter than minLen",
			seq:      strings.Repeat("A", 5),
			mRNAPlus: true,
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			g := geomFor(len(c.seq), c.mRNAPlus)
			got := findPolyA(c.seq, g, p, g.terminus)
			if got.found != c.wantFound {
				t.Fatalf("findPolyA(%q).found = %v, want %v", c.seq, got.found, c.wantFound)
			}
			if !c.wantFound {
				return
			}
			if got.start != c.wantStart || got.end != c.wantEnd || got.length != c.wantLen {
				t.Errorf("findPolyA(%q) = {start:%d end:%d len:%d}, want {start:%d end:%d len:%d}",
					c.seq, got.start, got.end, got.length, c.wantStart, c.wantEnd, c.wantLen)
			}
		})
	}
}

func TestPolyaRefPos(t *testing.T) {
	cases := []struct {
		name     string
		cigar    string
		pos      int
		startIdx int
		step     int
		wantPos  int
		wantSrc  string
		wantOk   bool
	}{
		{"trailing clip boundary", "50M20S", 100, 50, -1, 150, "clip", true},
		{"deep in trailing clip", "50M20S", 100, 60, -1, 150, "clip", true},
		{"aligned, sense", "50M20S", 100, 45, -1, 145, "aligned", true},
		{"leading clip boundary", "20S50M", 100, 19, 1, 99, "clip", true},
		{"aligned, antisense", "20S50M", 100, 25, 1, 105, "aligned", true},
		{"inside an insertion", "50M5I15S", 100, 52, -1, 149, "ins", true},
		{"past a deletion", "10M5D40M20S", 100, 30, -1, 135, "aligned", true},
		{"past a skip", "10M5N40M20S", 100, 30, -1, 135, "aligned", true},
		{"hard clip, trailing clip", "5H50M20S", 100, 50, -1, 150, "clip", true},
		{"hard clip, aligned", "5H50M20S", 100, 40, -1, 140, "aligned", true},
		{"no clips at all", "70M", 100, 60, -1, 160, "aligned", true},
		{"unspecified cigar", "*", 0, 5, -1, 0, "", false},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			seqLen := htsio.CigarQueryLen(c.cigar)
			if seqLen == 0 {
				seqLen = 10 // "*" cigar: SEQ length is irrelevant, but must be set
			}
			r := &htsio.SamRecord{
				Pos:   c.pos,
				Cigar: c.cigar,
				Seq:   strings.Repeat("A", seqLen),
			}
			pos, src, ok := polyaRefPos(r, c.startIdx, c.step)
			if ok != c.wantOk {
				t.Fatalf("polyaRefPos(%q@%d, %d) ok = %v, want %v", c.cigar, c.pos, c.startIdx, ok, c.wantOk)
			}
			if !ok {
				return
			}
			if pos != c.wantPos || src != c.wantSrc {
				t.Errorf("polyaRefPos(%q@%d, %d) = (%d, %q), want (%d, %q)",
					c.cigar, c.pos, c.startIdx, pos, src, c.wantPos, c.wantSrc)
			}
		})
	}
}

// polyaTestDefaults resets every ont-polya flag var to its registered default and
// restores it after the test, since they are package-level globals.
func polyaTestDefaults(t *testing.T) {
	t.Helper()
	reset := func() {
		polyaOutput = "-"
		polyaAntisense = false
		polyaMinAFrac = 0.8
		polyaMinLen = 10
		polyaWindow = 10
		polyaMaxJunk = 20
		polyaAdapter = ""
		polyaAdapterIdent = 0.75
		polyaShowLength = false
		polyaShowSource = false
		polyaTags = nil
		polyaNoNA = false
		polyaThreads = 1
		polyaCramRef = ""
	}
	reset()
	t.Cleanup(reset)
}

// runPolyaTSV writes records to a temp BAM, runs the command, and returns the
// output split into fields. Row 0 is the header.
func runPolyaTSV(t *testing.T, records []*htsio.SamRecord) [][]string {
	t.Helper()
	dir := t.TempDir()
	bamPath := filepath.Join(dir, "in.bam")
	outPath := filepath.Join(dir, "out.tsv")
	makeTestBAM(t, bamPath, records)

	polyaOutput = outPath
	if err := runPolya(bamPath); err != nil {
		t.Fatalf("runPolya: %v", err)
	}

	data, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	var rows [][]string
	for _, line := range strings.Split(strings.TrimRight(string(data), "\n"), "\n") {
		if line == "" {
			continue
		}
		rows = append(rows, strings.Split(line, "\t"))
	}
	return rows
}

// senseTail is 50 non-A bases followed by a 20 nt poly(A) tail.
func senseTail() string { return polyaNonA(50) + strings.Repeat("A", 20) }

// antisenseTail is a 20 nt tail stored as poly(T), followed by 50 non-T bases.
func antisenseTail() string { return strings.Repeat("T", 20) + polyaNonA(50) }

func TestPolya_StrandsAndSources(t *testing.T) {
	cases := []struct {
		name       string
		flag       int
		cigar      string
		pos        int
		seq        string
		wantPos    string
		wantStrand string
		wantSrc    string
	}{
		// Tail fully soft-clipped: the site is the first base past the alignment.
		{"sense clip", 0, "50M20S", 100, senseTail(), "150", "+", "clip"},
		// Tail absorbed into the alignment. pos 200 so a correct answer cannot
		// coincide with the sense-clip case's 150 and pass for the wrong reason.
		{"sense aligned", 0, "60M10S", 200, senseTail(), "250", "+", "aligned"},
		// Minus strand: the site is the base before the alignment.
		{"antisense clip", 16, "20S50M", 100, antisenseTail(), "99", "-", "clip"},
		{"antisense aligned", 16, "10S60M", 100, antisenseTail(), "109", "-", "aligned"},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			polyaTestDefaults(t)
			polyaShowLength = true
			polyaShowSource = true

			rows := runPolyaTSV(t, []*htsio.SamRecord{
				rec("r1", c.flag, "chr1", c.pos, c.cigar, c.seq, nil),
			})
			if len(rows) != 2 {
				t.Fatalf("got %d rows, want 2 (header + 1)", len(rows))
			}
			want := []string{"r1", "chr1", c.wantPos, c.wantStrand, "20", c.wantSrc}
			for i := range want {
				if rows[1][i] != want[i] {
					t.Errorf("row = %v, want %v", rows[1], want)
					break
				}
			}
		})
	}
}

func TestPolya_NoCall(t *testing.T) {
	polyaTestDefaults(t)

	rows := runPolyaTSV(t, []*htsio.SamRecord{
		rec("no_tail", 0, "chr1", 100, "70M", polyaNonA(70), nil),
	})
	if len(rows) != 2 {
		t.Fatalf("got %d rows, want 2", len(rows))
	}
	if rows[1][0] != "no_tail" || rows[1][2] != "NA" || rows[1][3] != "+" {
		t.Errorf("row = %v, want [no_tail chr1 NA +]", rows[1])
	}
}

func TestPolya_Unmapped(t *testing.T) {
	polyaTestDefaults(t)

	rows := runPolyaTSV(t, []*htsio.SamRecord{
		rec("unmapped", 4, "*", 0, "*", "*", nil),
	})
	if len(rows) != 2 {
		t.Fatalf("got %d rows, want 2", len(rows))
	}
	want := []string{"unmapped", "NA", "NA", "NA"}
	for i := range want {
		if rows[1][i] != want[i] {
			t.Errorf("row = %v, want %v", rows[1], want)
			break
		}
	}
}

func TestPolya_SecondarySupplementarySkipped(t *testing.T) {
	polyaTestDefaults(t)

	rows := runPolyaTSV(t, []*htsio.SamRecord{
		rec("sec", 0x100, "chr1", 100, "50M20S", senseTail(), nil),
		rec("supp", 0x800, "chr1", 100, "50M20S", senseTail(), nil),
	})
	if len(rows) != 1 {
		t.Fatalf("got %d rows, want 1 (header only): %v", len(rows), rows)
	}
}

func TestPolya_AntisenseOverride(t *testing.T) {
	polyaTestDefaults(t)
	polyaAntisense = true

	// FLAG 0 would normally mean a plus-strand mRNA, but --antisense inverts it,
	// so the tail is read as poly(T) in the leading clip.
	rows := runPolyaTSV(t, []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "20S50M", antisenseTail(), nil),
	})
	if len(rows) != 2 {
		t.Fatalf("got %d rows, want 2", len(rows))
	}
	if rows[1][2] != "99" || rows[1][3] != "-" {
		t.Errorf("row = %v, want pos 99 strand -", rows[1])
	}
}

func TestPolya_NoNA(t *testing.T) {
	polyaTestDefaults(t)

	records := []*htsio.SamRecord{
		rec("called", 0, "chr1", 100, "50M20S", senseTail(), nil),
		rec("no_tail", 0, "chr1", 200, "70M", polyaNonA(70), nil),
	}

	rows := runPolyaTSV(t, records)
	if len(rows) != 3 {
		t.Fatalf("default: got %d rows, want 3 (header + 2)", len(rows))
	}

	polyaNoNA = true
	rows = runPolyaTSV(t, records)
	if len(rows) != 2 {
		t.Fatalf("--no-na: got %d rows, want 2 (header + 1)", len(rows))
	}
	if rows[1][0] != "called" {
		t.Errorf("--no-na kept %q, want only the called read", rows[1][0])
	}
}

func TestPolya_TagColumns(t *testing.T) {
	polyaTestDefaults(t)
	polyaTags = []string{"pt"}

	rows := runPolyaTSV(t, []*htsio.SamRecord{
		rec("with_pt", 0, "chr1", 100, "50M20S", senseTail(), tags("pt", 'i', "123")),
		rec("no_pt", 0, "chr1", 200, "50M20S", senseTail(), nil),
	})
	if len(rows) != 3 {
		t.Fatalf("got %d rows, want 3", len(rows))
	}
	if rows[0][4] != "pt" {
		t.Errorf("header = %v, want a pt column at index 4", rows[0])
	}
	if rows[1][4] != "123" {
		t.Errorf("with_pt row = %v, want pt = 123", rows[1])
	}
	if rows[2][4] != "NA" {
		t.Errorf("no_pt row = %v, want pt = NA", rows[2])
	}
}

func TestPolya_OptionalColumns(t *testing.T) {
	polyaTestDefaults(t)

	records := []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "50M20S", senseTail(), nil),
	}

	rows := runPolyaTSV(t, records)
	want := []string{"read_name", "chrom", "polya_pos", "strand"}
	if strings.Join(rows[0], ",") != strings.Join(want, ",") {
		t.Errorf("default header = %v, want %v", rows[0], want)
	}

	polyaShowLength = true
	polyaShowSource = true
	rows = runPolyaTSV(t, records)
	want = []string{"read_name", "chrom", "polya_pos", "strand", "polya_len", "polya_source"}
	if strings.Join(rows[0], ",") != strings.Join(want, ",") {
		t.Errorf("header = %v, want %v", rows[0], want)
	}
}

func TestPolya_ContigClamp(t *testing.T) {
	polyaTestDefaults(t)

	// A minus-strand read at pos 1: the site would be pos-1 == 0, which is not a
	// valid 1-based coordinate.
	rows := runPolyaTSV(t, []*htsio.SamRecord{
		rec("r1", 16, "chr1", 1, "20S50M", antisenseTail(), nil),
	})
	if len(rows) != 2 {
		t.Fatalf("got %d rows, want 2", len(rows))
	}
	if rows[1][2] != "1" {
		t.Errorf("polya_pos = %q, want 1 (clamped from 0)", rows[1][2])
	}
}

func TestPolya_Adapter(t *testing.T) {
	polyaTestDefaults(t)
	polyaAdapter = "AGATCGGAAGAGC"
	polyaShowLength = true

	// An untrimmed adapter riding past the tail: 37 genomic + 20 nt tail + adapter.
	seq := polyaNonA(37) + strings.Repeat("A", 20) + "AGATCGGAAGAGC"
	rows := runPolyaTSV(t, []*htsio.SamRecord{
		rec("r1", 0, "chr1", 100, "37M33S", seq, nil),
	})
	if len(rows) != 2 {
		t.Fatalf("got %d rows, want 2", len(rows))
	}
	// The tail starts at SEQ index 37, in the soft clip, so the site is the first
	// base past the 37M alignment: 99 + 37 + 1.
	if rows[1][2] != "137" || rows[1][4] != "20" {
		t.Errorf("row = %v, want pos 137 and polya_len 20", rows[1])
	}
}
