package ontcmd

import (
	"testing"

	"github.com/compgen-io/cgltk/htsio"
)

func TestDetectSeparator(t *testing.T) {
	if sep := detectSeparator("AAAA-CCCC-GGGG-AAAA"); sep != "-" {
		t.Errorf("got %q, want %q", sep, "-")
	}
	if sep := detectSeparator("AAAATTCCCCTTGGGGTTAAAA"); sep != "TT" {
		t.Errorf("got %q, want %q", sep, "TT")
	}
}

func TestCountNonSepBases(t *testing.T) {
	tests := []struct {
		umi  string
		sep  string
		want int
	}{
		{"AAAA-CCCC-GGGG-AAAA", "-", 16},
		{"AAAATTCCCCTTGGGGTTAAAA", "TT", 16},
		{"AAAA-CCCC", "-", 8},
	}
	for _, tt := range tests {
		got := countNonSepBases(tt.umi, tt.sep)
		if got != tt.want {
			t.Errorf("countNonSepBases(%q, %q) = %d, want %d", tt.umi, tt.sep, got, tt.want)
		}
	}
}

func TestClusterUMIs_Identical(t *testing.T) {
	umiCounts := map[string]int{
		"AAAA-CCCC-GGGG-AAAA": 10,
	}
	canonical := make(map[string]string)
	clusterUMIs(umiCounts, canonical, false)

	if got := canonical["AAAA-CCCC-GGGG-AAAA"]; got != "AAAA-CCCC-GGGG-AAAA" {
		t.Errorf("canonical = %q, want self", got)
	}
}

func TestClusterUMIs_SimilarMerge(t *testing.T) {
	// Two UMIs that differ by 1 base (well within 12/16 threshold)
	umiCounts := map[string]int{
		"AAAA-CCCC-GGGG-AAAA": 10,
		"AAAA-CCCC-GGGG-AAAT": 3,
	}
	canonical := make(map[string]string)
	umiMergeMatchThreshold = 12
	clusterUMIs(umiCounts, canonical, false)

	if got := canonical["AAAA-CCCC-GGGG-AAAT"]; got != "AAAA-CCCC-GGGG-AAAA" {
		t.Errorf("canonical of variant = %q, want %q", got, "AAAA-CCCC-GGGG-AAAA")
	}
	if got := canonical["AAAA-CCCC-GGGG-AAAA"]; got != "AAAA-CCCC-GGGG-AAAA" {
		t.Errorf("canonical of most common = %q, want self", got)
	}
}

func TestClusterUMIs_DissimilarNoMerge(t *testing.T) {
	// Two UMIs that differ by many bases (should not merge at 12/16)
	umiCounts := map[string]int{
		"AAAA-AAAA-AAAA-AAAA": 10,
		"CCCC-CCCC-CCCC-CCCC": 5,
	}
	canonical := make(map[string]string)
	umiMergeMatchThreshold = 12
	clusterUMIs(umiCounts, canonical, false)

	if got := canonical["CCCC-CCCC-CCCC-CCCC"]; got != "CCCC-CCCC-CCCC-CCCC" {
		t.Errorf("dissimilar UMI should stay self, got %q", got)
	}
}

func TestClusterUMIs_TTSeparator(t *testing.T) {
	umiCounts := map[string]int{
		"AAAATTCCCCTTGGGGTTAAAA": 10,
		"AAAATTCCCCTTGGGGTTAAAT": 3,
	}
	canonical := make(map[string]string)
	umiMergeMatchThreshold = 12
	clusterUMIs(umiCounts, canonical, false)

	if got := canonical["AAAATTCCCCTTGGGGTTAAAT"]; got != "AAAATTCCCCTTGGGGTTAAAA" {
		t.Errorf("canonical of TT variant = %q, want %q", got, "AAAATTCCCCTTGGGGTTAAAA")
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
