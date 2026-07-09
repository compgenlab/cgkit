package ontcmd

import (
	"testing"

	"github.com/compgenlab/hts/htsio"
)

func TestDetectSeparator(t *testing.T) {
	if sep := detectSeparator("AAAA/CCCC/GGGG/AAAA"); sep != "/" {
		t.Errorf("got %q, want %q", sep, "/")
	}
	if sep := detectSeparator("AAAA-CCCC-GGGG-AAAA"); sep != "-" {
		t.Errorf("got %q, want %q", sep, "-")
	}
	if sep := detectSeparator("AAAATTCCCCTTGGGGTTAAAA"); sep != "TT" {
		t.Errorf("got %q, want %q", sep, "TT")
	}
}

func TestNormalizeUMISeparator(t *testing.T) {
	tests := []struct {
		umi  string
		want string
	}{
		{"AAAA/CCCC/GGGG/AAAA", "AAAA/CCCC/GGGG/AAAA"},
		{"AAAATTCCCCTTGGGGTTAAAA", "AAAA/CCCC/GGGG/AAAA"},
		{"AAAA-CCCC", "AAAA/CCCC"},
	}
	for _, tt := range tests {
		got := normalizeUMISeparator(tt.umi)
		if got != tt.want {
			t.Errorf("normalizeUMISeparator(%q) = %q, want %q", tt.umi, got, tt.want)
		}
	}
}

func TestUmiLevenshtein(t *testing.T) {
	tests := []struct {
		a, b string
		want int
	}{
		{"AAAA/CCCC/GGGG/AAAA", "AAAA/CCCC/GGGG/AAAA", 0},
		{"AAAA/CCCC/GGGG/AAAA", "AAAA/CCCC/GGGG/AAAT", 1},
		{"AAAA/CCCC/GGGG/AAAA", "AAAA/CCCC/GGGG/AACT", 2},
		// missed separator: edit distance 1 (insert "/")
		{"AAAA/CCCCGGGG/AAAA", "AAAA/CCCC/GGGG/AAAA", 1},
		// extra base in first group: edit distance 1
		{"AAAAA/CCCC/GGGG/AAAA", "AAAA/CCCC/GGGG/AAAA", 1},
		// completely different UMIs
		{"AAAA/AAAA/AAAA/AAAA", "CCCC/CCCC/CCCC/CCCC", 16},
	}
	for _, tt := range tests {
		var buf levBuf
		// Unbounded (-1) so the test exercises the full-distance path.
		got := levDist(tt.a, tt.b, &buf, -1)
		if got != tt.want {
			t.Errorf("levDist(%q, %q) = %d, want %d", tt.a, tt.b, got, tt.want)
		}
	}
}

func TestClusterUMIs_Identical(t *testing.T) {
	umiCounts := map[string]int{
		"AAAA/CCCC/GGGG/AAAA": 10,
	}
	representative := make(map[string]string)
	_, _ = clusterUMIs(umiCounts, representative, 1)

	if got := representative["AAAA/CCCC/GGGG/AAAA"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("representative = %q, want self", got)
	}
}

func TestClusterUMIs_SimilarMerge(t *testing.T) {
	umiCounts := map[string]int{
		"AAAA/CCCC/GGGG/AAAA": 10,
		"AAAA/CCCC/GGGG/AAAT": 3,
	}
	representative := make(map[string]string)
	umiClusterEditThreshold = 3
	results, _ := clusterUMIs(umiCounts, representative, 1)

	// Both cluster together; most common UMI (count=10) is picked as representative
	if got := representative["AAAA/CCCC/GGGG/AAAT"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("representative of variant = %q, want %q", got, "AAAA/CCCC/GGGG/AAAA")
	}
	if got := representative["AAAA/CCCC/GGGG/AAAA"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("representative of anchor = %q, want self", got)
	}

	// maxIntraClustDist should be 1 (one mismatch between the two UMIs)
	for _, r := range results {
		if r.maxIntraClustDist != 1 {
			t.Errorf("maxIntraClustDist = %d, want 1", r.maxIntraClustDist)
		}
	}
}

func TestClusterUMIs_DissimilarNoMerge(t *testing.T) {
	umiCounts := map[string]int{
		"AAAA/AAAA/AAAA/AAAA": 10,
		"CCCC/CCCC/CCCC/CCCC": 5,
	}
	representative := make(map[string]string)
	umiClusterEditThreshold = 3
	_, _ = clusterUMIs(umiCounts, representative, 1)

	if got := representative["CCCC/CCCC/CCCC/CCCC"]; got != "CCCC/CCCC/CCCC/CCCC" {
		t.Errorf("dissimilar UMI should stay self, got %q", got)
	}
}

func TestClusterUMIs_TTSeparator(t *testing.T) {
	umiCounts := map[string]int{
		"AAAATTCCCCTTGGGGTTAAAA": 10,
		"AAAATTCCCCTTGGGGTTAAAT": 3,
	}
	representative := make(map[string]string)
	umiClusterEditThreshold = 3
	_, _ = clusterUMIs(umiCounts, representative, 1)

	// Consensus is output in "/" form (normalized)
	if got := representative["AAAATTCCCCTTGGGGTTAAAT"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("representative of TT variant = %q, want %q", got, "AAAA/CCCC/GGGG/AAAA")
	}
}

func TestClusterUMIs_MissedSeparator(t *testing.T) {
	// xxxx-xxxxxxxx-xxxx (missed one separator) should cluster with xxxx-xxxx-xxxx-xxxx
	umiCounts := map[string]int{
		"AAAA/CCCC/GGGG/AAAA": 10,
		"AAAA-CCCCGGGG-AAAA":  3, // missing separator between CCCC and GGGG
	}
	representative := make(map[string]string)
	umiClusterEditThreshold = 3
	_, _ = clusterUMIs(umiCounts, representative, 1)

	if got := representative["AAAA-CCCCGGGG-AAAA"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("missed-separator UMI representative = %q, want %q", got, "AAAA/CCCC/GGGG/AAAA")
	}
}

func TestClusterUMIs_ExtraBase(t *testing.T) {
	// Most common UMI (count=10) is picked; extra-base variant clusters but doesn't win
	umiCounts := map[string]int{
		"AAAA/CCCC/GGGG/AAAA":  10,
		"GAAAA/CCCC/GGGG/AAAA": 3,
	}
	representative := make(map[string]string)
	umiClusterEditThreshold = 3
	_, _ = clusterUMIs(umiCounts, representative, 1)

	if got := representative["GAAAA/CCCC/GGGG/AAAA"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("extra-base UMI representative = %q, want %q", got, "AAAA/CCCC/GGGG/AAAA")
	}
}

func TestClusterUMIs_SingleLinkageChaining(t *testing.T) {
	// A and C are edit distance 6 apart (too far to cluster directly),
	// but B is within distance 3 of both. Single-linkage should merge all three.
	// A: AAAA/CCCC/GGGG/AAAA
	// B: AAAA/CCCC/GGGG/AAAT  (dist A-B = 1)
	// C: AAAA-CCCC-GGGG-TTTT  (dist B-C = 3, dist A-C = 4... let me pick carefully)
	//
	// A: AAAA/CCCC/GGGG/AAAA
	// B: AAAA-CCCC-GGGG-AATТ  (dist A-B = 2)
	// C: AAAA-CCCC-GGGG-TTTT  (dist B-C = 2, dist A-C = 4 > threshold=3)
	umiCounts := map[string]int{
		"AAAA/CCCC/GGGG/AAAA": 5,
		"AAAA-CCCC-GGGG-AATT": 5, // dist to A = 2, dist to C = 2
		"AAAA-CCCC-GGGG-TTTT": 5, // dist to A = 4 (> threshold), dist to B = 2
	}
	representative := make(map[string]string)
	umiClusterEditThreshold = 3
	_, _ = clusterUMIs(umiCounts, representative, 1)

	// All three should be in the same cluster (connected via B)
	ca := representative["AAAA/CCCC/GGGG/AAAA"]
	cb := representative["AAAA-CCCC-GGGG-AATT"]
	cc := representative["AAAA-CCCC-GGGG-TTTT"]
	if ca != cb || cb != cc {
		t.Errorf("single-linkage chaining failed: A=%q B=%q C=%q, all should be equal", ca, cb, cc)
	}
}

func TestClusterUMIs_AdaptiveMinDist(t *testing.T) {
	// Five short (4-base) UMIs. Only AAAA/AAAT are within edit distance
	// (d=1); the other three are mutually >3 apart, so the sole edge is
	// the AAAA-AAAT pair at d=1. With 10 total pairs and a 4-base UMI, the
	// expected random-collision rate at d=1 already blows past alpha=0.10,
	// so the adaptive filter *wants* to exclude d>=1 and collapse the
	// effective threshold to 0 (exact-match-only, no merge).
	umiCounts := map[string]int{
		"AAAA": 10,
		"AAAT": 3, // dist to AAAA = 1
		"CGCG": 5, // >3 from all others
		"GTGT": 5,
		"TCTC": 5,
	}

	saveEdit := umiClusterEditThreshold
	saveAdaptive := umiClusterAdaptiveThreshold
	saveAlpha := umiClusterAdaptiveAlpha
	saveMin := umiClusterAdaptiveMinDist
	t.Cleanup(func() {
		umiClusterEditThreshold = saveEdit
		umiClusterAdaptiveThreshold = saveAdaptive
		umiClusterAdaptiveAlpha = saveAlpha
		umiClusterAdaptiveMinDist = saveMin
	})

	umiClusterEditThreshold = 3
	umiClusterAdaptiveThreshold = true
	umiClusterAdaptiveAlpha = 0.10

	// Default floor (min=1): distance-1 edges are always kept, so AAAA and
	// AAAT still merge and the effective threshold floors at 1.
	umiClusterAdaptiveMinDist = 1
	rep := make(map[string]string)
	_, effThreshold := clusterUMIs(umiCounts, rep, 1)
	if effThreshold < 1 {
		t.Errorf("with --adaptive-min-dist 1, effective threshold = %d, want >= 1", effThreshold)
	}
	if got := rep["AAAT"]; got != "AAAA" {
		t.Errorf("with floor, AAAT representative = %q, want %q (should still merge)", got, "AAAA")
	}

	// Floor disabled (min=0): the filter is free to exclude d>=1, dropping
	// the only edge and preventing the merge.
	umiClusterAdaptiveMinDist = 0
	rep0 := make(map[string]string)
	_, effThreshold0 := clusterUMIs(umiCounts, rep0, 1)
	if effThreshold0 != 0 {
		t.Errorf("with --adaptive-min-dist 0, effective threshold = %d, want 0", effThreshold0)
	}
	if got := rep0["AAAT"]; got != "AAAT" {
		t.Errorf("without floor, AAAT representative = %q, want self (no merge)", got)
	}
}

func TestComputeRepresentativeUMI_MostCommon(t *testing.T) {
	// Highest count UMI is picked as representative
	members := []umiCount{
		{"AAAA/CCCC/GGGG/AAAA", 10},
		{"AAAA/CCCC/GGGG/AAAT", 3},
		{"AAAA/CCCC/GGGG/AAAT", 2},
	}
	got := computeRepresentativeUMI(members)
	if got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("computeRepresentativeUMI = %q, want %q", got, "AAAA/CCCC/GGGG/AAAA")
	}

	// When the other UMI has more reads, it wins
	members2 := []umiCount{
		{"AAAA/CCCC/GGGG/AAAA", 3},
		{"AAAA/CCCC/GGGG/AAAT", 8},
	}
	got2 := computeRepresentativeUMI(members2)
	if got2 != "AAAA/CCCC/GGGG/AAAT" {
		t.Errorf("computeRepresentativeUMI = %q, want %q", got2, "AAAA/CCCC/GGGG/AAAT")
	}

	// Tie broken by longer normalized length
	members3 := []umiCount{
		{"AAAA/CCCC/GGGG/AAAA", 5},
		{"GAAAA/CCCC/GGGG/AAAA", 5},
	}
	got3 := computeRepresentativeUMI(members3)
	if got3 != "GAAAA/CCCC/GGGG/AAAA" {
		t.Errorf("computeRepresentativeUMI tie = %q, want %q", got3, "GAAAA/CCCC/GGGG/AAAA")
	}
}

func TestUnionFind(t *testing.T) {
	uf := newUnionFind(10)

	// Initially all elements are their own root.
	for i := 0; i < 10; i++ {
		if got := uf.find(i); got != i {
			t.Errorf("find(%d) = %d, want %d", i, got, i)
		}
	}

	// Union 0 and 1.
	newRoot, oldRoot, merged := uf.union(0, 1)
	if !merged {
		t.Fatal("union(0,1) should merge")
	}
	if uf.find(0) != uf.find(1) {
		t.Error("0 and 1 should be in same set")
	}
	if newRoot == oldRoot {
		t.Error("newRoot and oldRoot should differ")
	}

	// Union 0 and 1 again — should be no-op.
	_, _, merged = uf.union(0, 1)
	if merged {
		t.Error("union(0,1) again should not merge")
	}

	// Transitive: union 1 and 2 should put 0,1,2 together.
	uf.union(1, 2)
	if uf.find(0) != uf.find(2) {
		t.Error("0 and 2 should be in same set after transitive union")
	}

	// 3 should still be separate.
	if uf.find(0) == uf.find(3) {
		t.Error("0 and 3 should be in different sets")
	}

	// Grow and use new elements.
	uf.grow(20)
	if got := uf.find(15); got != 15 {
		t.Errorf("find(15) after grow = %d, want 15", got)
	}
}

func TestEndProximityMatching_BothEnds(t *testing.T) {
	// AND mode: reads must match on BOTH 5' and 3' ends.
	gap := 50

	tests := []struct {
		name      string
		aStart    int
		aEnd      int
		bStart    int
		bEnd      int
		wantMatch bool
	}{
		{"identical positions", 100, 1000, 100, 1000, true},
		{"5' within gap, 3' within gap", 100, 1000, 140, 1030, true},
		{"5' within gap, 3' outside gap", 100, 1000, 140, 1100, false},
		{"5' at exact gap boundary, 3' match", 100, 1000, 150, 1000, true},
		{"5' beyond gap, 3' match", 100, 1000, 151, 1000, false},
		{"same start, different length (3' far)", 100, 1000, 100, 500, false},
		{"truncated read, 3' close", 100, 1000, 120, 980, true},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Simulate AND mode matching: 5' within gap AND 3' within gap
			fivePrime := (tt.bStart-tt.aStart) <= gap && (tt.aStart-tt.bStart) <= gap
			diff := tt.bEnd - tt.aEnd
			if diff < 0 {
				diff = -diff
			}
			threePrime := diff <= gap
			got := fivePrime && threePrime
			if got != tt.wantMatch {
				t.Errorf("AND match = %v, want %v (5'=%v, 3'=%v)", got, tt.wantMatch, fivePrime, threePrime)
			}
		})
	}
}

func TestEndProximityMatching_SingleEnd(t *testing.T) {
	// OR mode: reads match if EITHER 5' or 3' ends are within gap.
	gap := 50

	tests := []struct {
		name      string
		aStart    int
		aEnd      int
		bStart    int
		bEnd      int
		wantMatch bool
	}{
		{"identical positions", 100, 1000, 100, 1000, true},
		{"5' within gap, 3' outside gap", 100, 1000, 140, 1100, true},
		{"5' outside gap, 3' within gap", 100, 1000, 200, 1030, true},
		{"both outside gap", 100, 1000, 200, 1100, false},
		{"5' at boundary, 3' far", 100, 1000, 150, 2000, true},
		{"truncated molecule (5' match, 3' far)", 100, 1000, 120, 500, true},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fivePrime := (tt.bStart-tt.aStart) <= gap && (tt.aStart-tt.bStart) <= gap
			diff := tt.bEnd - tt.aEnd
			if diff < 0 {
				diff = -diff
			}
			threePrime := diff <= gap
			got := fivePrime || threePrime
			if got != tt.wantMatch {
				t.Errorf("OR match = %v, want %v (5'=%v, 3'=%v)", got, tt.wantMatch, fivePrime, threePrime)
			}
		})
	}
}

func TestEjectionSafety_BothEnds(t *testing.T) {
	// AND mode: eject when curStart - B.start > gap
	gap := 50

	tests := []struct {
		name        string
		bStart      int
		curStart    int
		wantEjected bool
	}{
		{"within gap", 100, 140, false},
		{"at gap boundary", 100, 150, false},
		{"beyond gap", 100, 151, true},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.curStart-tt.bStart > gap
			if got != tt.wantEjected {
				t.Errorf("ejected = %v, want %v", got, tt.wantEjected)
			}
		})
	}
}

func TestEjectionSafety_SingleEnd(t *testing.T) {
	// OR mode: eject when curStart > B.end + gap
	gap := 50

	tests := []struct {
		name        string
		bEnd        int
		curStart    int
		wantEjected bool
	}{
		{"well before end", 1000, 500, false},
		{"at end+gap boundary", 1000, 1050, false},
		{"beyond end+gap", 1000, 1051, true},
		{"long read stays buffered", 500000, 1000, false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.curStart > tt.bEnd+gap
			if got != tt.wantEjected {
				t.Errorf("ejected = %v, want %v", got, tt.wantEjected)
			}
		})
	}
}

func TestValidateCoordinateSorted(t *testing.T) {
	tests := []struct {
		name    string
		lines   []string
		wantErr bool
	}{
		{
			"coordinate sorted",
			[]string{"@HD\tVN:1.6\tSO:coordinate"},
			false,
		},
		{
			"unsorted",
			[]string{"@HD\tVN:1.6\tSO:unsorted"},
			true,
		},
		{
			"no HD line",
			[]string{"@SQ\tSN:chr1\tLN:100"},
			true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			h := htsio.NewSamHeader()
			for _, line := range tt.lines {
				h.AddLine(line)
			}
			err := validateCoordinateSorted(h)
			if (err != nil) != tt.wantErr {
				t.Errorf("validateCoordinateSorted() err = %v, wantErr = %v", err, tt.wantErr)
			}
		})
	}
}

func TestExtractJunctions(t *testing.T) {
	tests := []struct {
		name     string
		cigar    string
		refStart int
		want     []spliceJunction
	}{
		{"no junctions", "100M", 0, nil},
		{"single junction", "50M1000N50M", 100, []spliceJunction{{150, 1150}}},
		{"two junctions", "30M500N20M200N40M", 0, []spliceJunction{{30, 530}, {550, 750}}},
		{"junction with soft clip", "5S50M1000N50M", 100, []spliceJunction{{150, 1150}}},
		{"junction with insertion", "30M5I20M1000N50M", 0, []spliceJunction{{50, 1050}}},
		{"junction with deletion", "30M5D20M1000N50M", 0, []spliceJunction{{55, 1055}}},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := extractJunctions(tt.cigar, tt.refStart)
			if len(got) != len(tt.want) {
				t.Fatalf("extractJunctions(%q, %d) = %v, want %v", tt.cigar, tt.refStart, got, tt.want)
			}
			for i := range got {
				if got[i] != tt.want[i] {
					t.Errorf("junction[%d] = %v, want %v", i, got[i], tt.want[i])
				}
			}
		})
	}
}

func TestMergeAdjacentJunctions(t *testing.T) {
	tests := []struct {
		name   string
		input  []spliceJunction
		window int
		want   []spliceJunction
	}{
		{"nil", nil, 10, nil},
		{"single", []spliceJunction{{100, 200}}, 10, []spliceJunction{{100, 200}}},
		{"far apart", []spliceJunction{{100, 200}, {300, 400}}, 10, []spliceJunction{{100, 200}, {300, 400}}},
		{"small exon merged", []spliceJunction{{100, 200}, {208, 400}}, 10, []spliceJunction{{100, 400}}},
		{"exact window", []spliceJunction{{100, 200}, {210, 400}}, 10, []spliceJunction{{100, 400}}},
		{"just outside window", []spliceJunction{{100, 200}, {211, 400}}, 10, []spliceJunction{{100, 200}, {211, 400}}},
		{"three merged to one", []spliceJunction{{100, 200}, {205, 300}, {305, 500}}, 10, []spliceJunction{{100, 500}}},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := mergeAdjacentJunctions(tt.input, tt.window)
			if len(got) != len(tt.want) {
				t.Fatalf("mergeAdjacentJunctions() = %v, want %v", got, tt.want)
			}
			for i := range got {
				if got[i] != tt.want[i] {
					t.Errorf("junction[%d] = %v, want %v", i, got[i], tt.want[i])
				}
			}
		})
	}
}

func TestJunctionsCompatible(t *testing.T) {
	mkRead := func(start, end int, junctions []spliceJunction) *bufferedRead {
		return &bufferedRead{start: start, end: end, junctions: junctions}
	}

	t.Run("both empty", func(t *testing.T) {
		a := mkRead(0, 100, nil)
		b := mkRead(0, 100, nil)
		if !junctionsCompatible(a, b, 10, false, 50) {
			t.Error("both empty should be compatible")
		}
	})

	t.Run("one empty one not", func(t *testing.T) {
		a := mkRead(0, 100, nil)
		b := mkRead(0, 100, []spliceJunction{{50, 150}})
		if junctionsCompatible(a, b, 10, false, 50) {
			t.Error("one empty, one not should be incompatible")
		}
	})

	t.Run("exact match default mode", func(t *testing.T) {
		a := mkRead(0, 200, []spliceJunction{{50, 150}})
		b := mkRead(0, 200, []spliceJunction{{55, 148}})
		if !junctionsCompatible(a, b, 10, false, 50) {
			t.Error("within window should match")
		}
	})

	t.Run("outside window default mode", func(t *testing.T) {
		a := mkRead(0, 200, []spliceJunction{{50, 150}})
		b := mkRead(0, 200, []spliceJunction{{65, 150}})
		if junctionsCompatible(a, b, 10, false, 50) {
			t.Error("outside window should not match")
		}
	})

	t.Run("different count default mode", func(t *testing.T) {
		a := mkRead(0, 200, []spliceJunction{{50, 150}})
		b := mkRead(0, 200, []spliceJunction{{50, 150}, {200, 300}})
		if junctionsCompatible(a, b, 10, false, 50) {
			t.Error("different junction count should not match in default mode")
		}
	})

	t.Run("match-one-end suffix (3' anchor)", func(t *testing.T) {
		// Reads share 3' end, shorter read has suffix of junctions
		a := mkRead(0, 1000, []spliceJunction{{100, 200}, {400, 600}, {700, 900}})
		b := mkRead(300, 1000, []spliceJunction{{400, 600}, {700, 900}})
		if !junctionsCompatible(a, b, 10, true, 50) {
			t.Error("suffix match with 3' anchor should be compatible")
		}
	})

	t.Run("match-one-end prefix (5' anchor)", func(t *testing.T) {
		// Reads share 5' start, shorter read has prefix of junctions
		a := mkRead(0, 1000, []spliceJunction{{100, 200}, {400, 600}, {700, 900}})
		b := mkRead(0, 650, []spliceJunction{{100, 200}, {400, 600}})
		if !junctionsCompatible(a, b, 10, true, 50) {
			t.Error("prefix match with 5' anchor should be compatible")
		}
	})

	t.Run("match-one-end no anchor match", func(t *testing.T) {
		// Neither end matches within overlap
		a := mkRead(0, 1000, []spliceJunction{{100, 200}, {400, 600}, {700, 900}})
		b := mkRead(200, 800, []spliceJunction{{400, 600}})
		if junctionsCompatible(a, b, 10, true, 50) {
			t.Error("no anchor match should be incompatible")
		}
	})
}
