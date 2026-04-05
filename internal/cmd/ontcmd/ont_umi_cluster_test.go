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
		got := umiLevenshtein(tt.a, tt.b)
		if got != tt.want {
			t.Errorf("umiLevenshtein(%q, %q) = %d, want %d", tt.a, tt.b, got, tt.want)
		}
	}
}

func TestClusterUMIs_Identical(t *testing.T) {
	umiCounts := map[string]int{
		"AAAA/CCCC/GGGG/AAAA": 10,
	}
	consensus := make(map[string]string)
	clusterUMIs(umiCounts, consensus, false)

	if got := consensus["AAAA/CCCC/GGGG/AAAA"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("consensus = %q, want self", got)
	}
}

func TestClusterUMIs_SimilarMerge(t *testing.T) {
	umiCounts := map[string]int{
		"AAAA/CCCC/GGGG/AAAA": 10,
		"AAAA/CCCC/GGGG/AAAT": 3,
	}
	consensus := make(map[string]string)
	umiClusterEditThreshold = 3
	results := clusterUMIs(umiCounts, consensus, false)

	// Both cluster together; most common UMI (count=10) is picked as representative
	if got := consensus["AAAA/CCCC/GGGG/AAAT"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("consensus of variant = %q, want %q", got, "AAAA/CCCC/GGGG/AAAA")
	}
	if got := consensus["AAAA/CCCC/GGGG/AAAA"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("consensus of anchor = %q, want self", got)
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
	consensus := make(map[string]string)
	umiClusterEditThreshold = 3
	clusterUMIs(umiCounts, consensus, false)

	if got := consensus["CCCC/CCCC/CCCC/CCCC"]; got != "CCCC/CCCC/CCCC/CCCC" {
		t.Errorf("dissimilar UMI should stay self, got %q", got)
	}
}

func TestClusterUMIs_TTSeparator(t *testing.T) {
	umiCounts := map[string]int{
		"AAAATTCCCCTTGGGGTTAAAA": 10,
		"AAAATTCCCCTTGGGGTTAAAT": 3,
	}
	consensus := make(map[string]string)
	umiClusterEditThreshold = 3
	clusterUMIs(umiCounts, consensus, false)

	// Consensus is output in "/" form (normalized)
	if got := consensus["AAAATTCCCCTTGGGGTTAAAT"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("consensus of TT variant = %q, want %q", got, "AAAA/CCCC/GGGG/AAAA")
	}
}

func TestClusterUMIs_MissedSeparator(t *testing.T) {
	// xxxx-xxxxxxxx-xxxx (missed one separator) should cluster with xxxx-xxxx-xxxx-xxxx
	umiCounts := map[string]int{
		"AAAA/CCCC/GGGG/AAAA": 10,
		"AAAA-CCCCGGGG-AAAA":  3, // missing separator between CCCC and GGGG
	}
	consensus := make(map[string]string)
	umiClusterEditThreshold = 3
	clusterUMIs(umiCounts, consensus, false)

	if got := consensus["AAAA-CCCCGGGG-AAAA"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("missed-separator UMI consensus = %q, want %q", got, "AAAA/CCCC/GGGG/AAAA")
	}
}

func TestClusterUMIs_ExtraBase(t *testing.T) {
	// Most common UMI (count=10) is picked; extra-base variant clusters but doesn't win
	umiCounts := map[string]int{
		"AAAA/CCCC/GGGG/AAAA":  10,
		"GAAAA/CCCC/GGGG/AAAA": 3,
	}
	consensus := make(map[string]string)
	umiClusterEditThreshold = 3
	clusterUMIs(umiCounts, consensus, false)

	if got := consensus["GAAAA/CCCC/GGGG/AAAA"]; got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("extra-base UMI consensus = %q, want %q", got, "AAAA/CCCC/GGGG/AAAA")
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
	consensus := make(map[string]string)
	umiClusterEditThreshold = 3
	clusterUMIs(umiCounts, consensus, false)

	// All three should be in the same cluster (connected via B)
	ca := consensus["AAAA/CCCC/GGGG/AAAA"]
	cb := consensus["AAAA-CCCC-GGGG-AATT"]
	cc := consensus["AAAA-CCCC-GGGG-TTTT"]
	if ca != cb || cb != cc {
		t.Errorf("single-linkage chaining failed: A=%q B=%q C=%q, all should be equal", ca, cb, cc)
	}
}

func TestComputeConsensusUMI_MostCommon(t *testing.T) {
	// Highest count UMI is picked as representative
	members := []umiCount{
		{"AAAA/CCCC/GGGG/AAAA", 10},
		{"AAAA/CCCC/GGGG/AAAT", 3},
		{"AAAA/CCCC/GGGG/AAAT", 2},
	}
	got := computeConsensusUMI(members)
	if got != "AAAA/CCCC/GGGG/AAAA" {
		t.Errorf("computeConsensusUMI = %q, want %q", got, "AAAA/CCCC/GGGG/AAAA")
	}

	// When the other UMI has more reads, it wins
	members2 := []umiCount{
		{"AAAA/CCCC/GGGG/AAAA", 3},
		{"AAAA/CCCC/GGGG/AAAT", 8},
	}
	got2 := computeConsensusUMI(members2)
	if got2 != "AAAA/CCCC/GGGG/AAAT" {
		t.Errorf("computeConsensusUMI = %q, want %q", got2, "AAAA/CCCC/GGGG/AAAT")
	}

	// Tie broken by longer normalized length
	members3 := []umiCount{
		{"AAAA/CCCC/GGGG/AAAA", 5},
		{"GAAAA/CCCC/GGGG/AAAA", 5},
	}
	got3 := computeConsensusUMI(members3)
	if got3 != "GAAAA/CCCC/GGGG/AAAA" {
		t.Errorf("computeConsensusUMI tie = %q, want %q", got3, "GAAAA/CCCC/GGGG/AAAA")
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
