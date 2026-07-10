package vcfcmd

import (
	"strings"
	"testing"
)

func TestCanonicalGT(t *testing.T) {
	cases := map[string]string{
		"0/0": "0/0",
		"0/1": "0/1",
		"1/0": "0/1",
		"0|1": "0/1",
		"1|0": "0/1",
		"1/1": "1/1",
		"2/1": "1/2",
		"0/2": "0/2",
		"./.": "./.",
		".|.": "./.",
		"1/.": "./1",
		".|1": "./1",
		"1":   "1",
		".":   ".",
		"":    ".",
	}
	for raw, want := range cases {
		if got := canonicalGT(raw); got != want {
			t.Errorf("canonicalGT(%q) = %q, want %q", raw, got, want)
		}
	}
}

func TestGtColumnOrder(t *testing.T) {
	cols := []string{"./.", "1/2", "0/0", "./1", "1/1", "0/1", "0/2"}
	want := []string{"0/0", "0/1", "1/1", "0/2", "1/2", "./1", "./."}
	sortGtColumns(cols)
	if strings.Join(cols, ",") != strings.Join(want, ",") {
		t.Errorf("column order = %v, want %v", cols, want)
	}
}

func TestParseLocus(t *testing.T) {
	s, err := parseLocus("chr1:100")
	if err != nil || s.chrom != "chr1" || s.pos != 100 || s.hasRA {
		t.Fatalf("parseLocus(chr1:100) = %+v, %v", s, err)
	}
	s, err = parseLocus("chr1:300:G:GA")
	if err != nil || s.ref != "G" || s.alt != "GA" || !s.hasRA {
		t.Fatalf("parseLocus(chr1:300:G:GA) = %+v, %v", s, err)
	}
	if _, err := parseLocus("chr1"); err == nil {
		t.Error("parseLocus(chr1) expected error")
	}
	if _, err := parseLocus("chr1:abc"); err == nil {
		t.Error("parseLocus(chr1:abc) expected error")
	}
}

// dataLines returns the non-comment (non "##") output lines of vcf-gtcount.
func dataLines(out string) []string {
	var lines []string
	for _, l := range strings.Split(strings.TrimRight(out, "\n"), "\n") {
		if l != "" && !strings.HasPrefix(l, "##") {
			lines = append(lines, l)
		}
	}
	return lines
}

func TestVcfGtCount(t *testing.T) {
	out := runVcf(t, "vcf-gtcount", "testdata/sample.vcf.gz",
		"chr1:100", "chr1:300:G:GA", "chr9:999")
	lines := dataLines(out)
	want := []string{
		"chrom\tpos\tref\talt\tgt\tcount",
		"chr1\t100\tA\tG\t0/0\t1",
		"chr1\t100\tA\tG\t0/1\t1",
		"chr1\t300\tG\tGA\t0/0\t1",
		"chr1\t300\tG\tGA\t0/1\t1",
		"chr9\t999\t.\t.\t.\t0", // absent contig -> sentinel row
	}
	if strings.Join(lines, "\n") != strings.Join(want, "\n") {
		t.Errorf("output:\n%s\nwant:\n%s", strings.Join(lines, "\n"), strings.Join(want, "\n"))
	}
}

func TestVcfGtCountRefAltMismatch(t *testing.T) {
	// chr1:100 is A>G; requiring A>C must not match -> sentinel row.
	out := runVcf(t, "vcf-gtcount", "testdata/sample.vcf.gz", "chr1:100:A:C")
	lines := dataLines(out)
	if len(lines) != 2 || lines[1] != "chr1\t100\tA\tC\t.\t0" {
		t.Errorf("ref/alt mismatch output:\n%s", strings.Join(lines, "\n"))
	}
}

func TestVcfGtCountPassing(t *testing.T) {
	// chr1:200 is FILTER=lowqual; --passing drops it -> sentinel row.
	out := runVcf(t, "vcf-gtcount", "--passing", "testdata/sample.vcf.gz", "chr1:200")
	lines := dataLines(out)
	if len(lines) != 2 || lines[1] != "chr1\t200\t.\t.\t.\t0" {
		t.Errorf("--passing output:\n%s", strings.Join(lines, "\n"))
	}
}

func TestVcfGtCountThreadsOrder(t *testing.T) {
	// Several loci (repeated to give the pool work); -t 4 must produce the same
	// rows in the same site order as the serial default.
	loci := []string{
		"chr1:100", "chr1:200", "chr1:300", "chr2:500", "chr2:1000",
		"chr1:100", "chr9:999", "chr1:300", "chr1:200",
	}
	serial := dataLines(runVcf(t, append([]string{"vcf-gtcount", "testdata/sample.vcf.gz"}, loci...)...))
	parallel := dataLines(runVcf(t, append([]string{"vcf-gtcount", "-t", "4", "testdata/sample.vcf.gz"}, loci...)...))
	if strings.Join(serial, "\n") != strings.Join(parallel, "\n") {
		t.Errorf("-t 4 output differs from serial:\nserial:\n%s\nparallel:\n%s",
			strings.Join(serial, "\n"), strings.Join(parallel, "\n"))
	}
}

func TestGtProgressLine(t *testing.T) {
	cases := []struct {
		done, total int64
		want        string
	}{
		{0, 20, "vcf-gtcount: 0/20 sites (0%)"},
		{5, 20, "vcf-gtcount: 5/20 sites (25%)"},
		{20, 20, "vcf-gtcount: 20/20 sites (100%)"},
		{1, 0, "vcf-gtcount: 1/0 sites (0%)"}, // guard against divide-by-zero
	}
	for _, c := range cases {
		if got := gtProgressLine(c.done, c.total); got != c.want {
			t.Errorf("gtProgressLine(%d,%d) = %q, want %q", c.done, c.total, got, c.want)
		}
	}
}

func TestVcfGtCountThreadsInvalid(t *testing.T) {
	for _, n := range []string{"0", "-1"} {
		if err := runVcfErr(t, "vcf-gtcount", "-t", n, "testdata/sample.vcf.gz", "chr1:100"); err == nil {
			t.Errorf("-t %s expected an error", n)
		}
	}
}
