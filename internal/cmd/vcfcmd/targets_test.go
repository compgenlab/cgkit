package vcfcmd

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// Tests for --variant, the single site selector: inline loci, inline regions, or a
// file whose format is detected from its content.

// writeFile drops a fixture into a temp dir and returns its path.
func writeFile(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

// positions returns the pos column of every data row.
func positions(out string) []string {
	var got []string
	for _, r := range tsvDataRows(dataRowsOnly(out)) {
		got = append(got, strings.Split(r, "\t")[colPos])
	}
	return got
}

// TestTargetFormatDetection pins the rule that decides what a target file is.
//
// A VCF announces itself. Otherwise column 3 decides, because a BED's third
// column is an end coordinate and always numeric where a site list's is a REF
// allele and never is.
func TestTargetFormatDetection(t *testing.T) {
	cases := []struct {
		name    string
		content string
		want    targetFormat
	}{
		{"bed.txt", "chr1\t99\t300\n", targetBED},
		{"bed-extra-cols.txt", "chr1 99 300 name 0 +\n", targetBED},
		{"bed-with-track.txt", "track name=x\nchr1 99 300\n", targetBED},
		{"sites.txt", "chr1 100 A G\n", targetList},
		{"sites-extra-cols.txt", "chr1 100 A G rs123\n", targetList},
		{"sites-pos-only.txt", "chr1 100\n", targetList},
		{"sites-colon.txt", "chr1:100:A:G\n", targetList},
		{"sites-commented.txt", "# a comment\n\nchr1 100 A G\n", targetList},
		{"panel.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\n", targetVCF},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := detectTargetFile(writeFile(t, tc.name, tc.content))
			if err != nil {
				t.Fatalf("detect: %v", err)
			}
			if got != tc.want {
				t.Errorf("detected %q, want %q", got, tc.want)
			}
		})
	}
}

// TestTargetBedIsZeroBased is the coordinate check, and the reason detection has
// to be right rather than merely plausible: a BED read as a site list does not
// fail, it shifts every position by one.
func TestTargetBedIsZeroBased(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	// 0-based half-open [99,300) is 1-based 100..300.
	bed := writeFile(t, "p.bed", "chr1\t99\t300\n")
	got := positions(runVcf(t, "vcf-varquery", "--variant", bed,
		"--sample", "S1", "--hom-ref", "--min-dp", "10", base))
	want := "100,200,300"
	if strings.Join(got, ",") != want {
		t.Errorf("BED [99,300) selected %v, want %s -- if this is 99 or omits 100 the "+
			"file was read as 1-based", got, want)
	}

	// The same interval written as a site list means one position, not a range.
	list := writeFile(t, "p.txt", "chr1 99 300\n")
	if f, err := detectTargetFile(list); err != nil || f != targetBED {
		t.Fatalf("a three-column numeric line is a BED by construction; got %q %v", f, err)
	}
}

// TestTargetInlineForms pins the command-line grammar.
func TestTargetInlineForms(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	cases := []struct{ arg, want string }{
		{"chr1:500:A:T", "500"},     // exact locus
		{"chr1:300", "300"},         // any variant at a position
		{"chr1:250-450", "300,400"}, // a region
		{"chr2", "100"},             // a whole contig
	}
	for _, tc := range cases {
		got := positions(runVcf(t, "vcf-varquery", "--variant", tc.arg,
			"--hom-ref", "--min-dp", "10", base))
		if strings.Join(uniq(got), ",") != tc.want {
			t.Errorf("--variant %s selected %v, want %s", tc.arg, uniq(got), tc.want)
		}
	}
}

// TestTargetInlineErrorsAreAboutLoci pins the reason os.Stat is consulted first: a
// mistyped locus is not a file, so it must not be reported as a missing file.
func TestTargetInlineErrorsAreAboutLoci(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	err := runVcfErr(t, "vcf-varquery", "--variant", "chr1:100:A", base)
	if err == nil {
		t.Fatal("expected an error for a three-field locus")
	}
	if strings.Contains(err.Error(), "no such file") {
		t.Errorf("a mistyped locus should not be reported as a missing file: %v", err)
	}
	if !strings.Contains(err.Error(), "chr1:100:A") {
		t.Errorf("the error should name the bad value, got %v", err)
	}
}

// TestTargetFileAndInlineCombine pins that files and inline selectors accumulate
// rather than one replacing the other.
func TestTargetFileAndInlineCombine(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	list := writeFile(t, "p.txt", "chr1 100 A G\n")
	got := uniq(positions(runVcf(t, "vcf-varquery",
		"--variant", list, "--variant", "chr1:500:A:T", "--hom-ref", "--min-dp", "10", base)))
	if strings.Join(got, ",") != "100,500" {
		t.Errorf("a file plus an inline locus should give both: %v", got)
	}
}

// TestTargetEmptyFileIsRefused pins that an undecidable file errors rather than
// silently selecting nothing, which would read as a real negative result.
func TestTargetEmptyFileIsRefused(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	for _, content := range []string{"", "\n\n", "# only comments\n"} {
		err := runVcfErr(t, "vcf-varquery", "--variant",
			writeFile(t, "empty.txt", content), base)
		if err == nil {
			t.Errorf("an empty target file (%q) should be refused", content)
		}
	}
}

// TestSiteListWorksForBothCommands is the protection against format drift. The
// target parser here and vcf-gtcount's readSitesFile are separate code; what must
// not diverge is the file FORMAT, so one file has to parse the same through both.
//
// Compared at the parser level rather than by running both commands, because
// vcf-gtcount needs a tabix index for random access and that is irrelevant to the
// claim being made.
func TestSiteListWorksForBothCommands(t *testing.T) {
	list := writeFile(t, "shared.txt",
		"# shared site list\nchr1 100 A G\nchr1 500 A T\nchr2\t300\tG\tC\n")

	mine, err := parseTargets([]string{list})
	if err != nil {
		t.Fatalf("varquery parser: %v", err)
	}
	theirs, err := readSitesFile(list)
	if err != nil {
		t.Fatalf("vcf-gtcount parser: %v", err)
	}
	if len(mine.Loci) != len(theirs) {
		t.Fatalf("varquery read %d loci, vcf-gtcount %d, from the same file",
			len(mine.Loci), len(theirs))
	}
	for i, l := range mine.Loci {
		g := theirs[i]
		if l.Chrom != g.chrom || int(l.Pos) != g.pos || l.Ref != g.ref || l.Alt != g.alt {
			t.Errorf("row %d differs: varquery %s:%d:%s:%s, gtcount %s:%d:%s:%s",
				i, l.Chrom, l.Pos, l.Ref, l.Alt, g.chrom, g.pos, g.ref, g.alt)
		}
	}
}

// TestRegionFlagIsGone pins that --region was removed rather than left as a
// second way to say the same thing.
func TestRegionFlagIsGone(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	if err := runVcfErr(t, "vcf-varquery", "--region", "chr1", base); err == nil {
		t.Error("--region should no longer exist; regions are --variant chr1:start-end")
	}
}

func uniq(in []string) []string {
	seen := map[string]bool{}
	var out []string
	for _, s := range in {
		if !seen[s] {
			seen[s] = true
			out = append(out, s)
		}
	}
	return out
}
