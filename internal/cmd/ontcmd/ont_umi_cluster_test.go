package ontcmd

import (
	"testing"

	"github.com/compgen-io/cgltk/htsio"
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
			fivePrime := (tt.bStart - tt.aStart) <= gap && (tt.aStart - tt.bStart) <= gap
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
			fivePrime := (tt.bStart - tt.aStart) <= gap && (tt.aStart - tt.bStart) <= gap
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
