package vcfcmd

import (
	"context"
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
			got, err := detectTargetFile(context.Background(), writeFile(t, tc.name, tc.content))
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
	if f, err := detectTargetFile(context.Background(), list); err != nil || f != targetBED {
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

// TestTargetColonContigsAreNotLoci pins the disambiguation that lets a contig
// name carry colons.
//
// GRCh38's ALT contigs are named like HLA-A*01:01:01:01, which splits into four
// colon-fields exactly as a locus does. What separates them is that a locus's last
// two fields are REF and ALT alleles and never numeric, where the contig's are.
// Without that test the contig parses as chrom=HLA-A*01, pos=1, ref=01, alt=01:
// accepted, wrong, and silent.
func TestTargetColonContigsAreNotLoci(t *testing.T) {
	for _, name := range []string{
		"HLA-A*01:01:01:01", // four fields, all-numeric tail
		"HLA-B*07:02:01",    // three fields
		"chr1:100:A",        // a locus missing a field lands here too
	} {
		set, err := parseTargets(context.Background(), []string{name})
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if len(set.Loci) != 0 {
			t.Errorf("%s parsed as a locus %v; it is a contig name", name, set.Loci)
		}
		if len(set.Spans) != 1 || set.Spans[0].Chrom != name {
			t.Errorf("%s should be one whole-contig span, got %+v", name, set.Spans)
		}
	}

	// A real locus is still a locus, including an indel whose alleles are long.
	set, err := parseTargets(context.Background(), []string{"chr1:300:G:GATTACA"})
	if err != nil {
		t.Fatal(err)
	}
	if len(set.Loci) != 1 || set.Loci[0].Alt != "GATTACA" {
		t.Errorf("want one locus with alt GATTACA, got %+v / %+v", set.Loci, set.Spans)
	}
}

// TestTargetUnmatchedIsReported pins the safety net for that fallback: because a
// malformed locus becomes a contig nobody has, selecting nothing has to be said
// out loud rather than read as a real negative result.
func TestTargetUnmatchedIsReported(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	out := runVcf(t, "vcf-varquery", "--variant", "chr1:100:A", base)
	if !strings.Contains(out, "no rows for any target") {
		t.Errorf("a target that matched nothing should be reported:\n%s", out)
	}
	if !strings.Contains(out, "read as a contig name") {
		t.Errorf("the warning should say how the value was interpreted:\n%s", out)
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

// TestSiteListWorksForBothCommands pins that one file serves both commands.
//
// They now share a single parser rather than two that had to agree, so this checks
// the remaining seam: vcf-gtcount's conversion of the shared target set into its
// own site type.
func TestSiteListWorksForBothCommands(t *testing.T) {
	list := writeFile(t, "shared.txt",
		"# shared site list\nchr1 100 A G\nchr1 500 A T\nchr2:300:G:C\n")

	targets, err := parseTargets(context.Background(), []string{list})
	if err != nil {
		t.Fatalf("target parser: %v", err)
	}
	if len(targets.Loci) != 3 {
		t.Fatalf("want 3 loci from the shared list, got %d", len(targets.Loci))
	}

	sites, err := collectGtSites(context.Background(), nil, []string{list})
	if err != nil {
		t.Fatalf("vcf-gtcount rejected the shared list: %v", err)
	}
	if len(sites) != len(targets.Loci) {
		t.Fatalf("varquery read %d targets, vcf-gtcount %d", len(targets.Loci), len(sites))
	}
	for i, l := range targets.Loci {
		g := sites[i]
		if l.Chrom != g.chrom || int(l.Pos) != g.pos || l.Ref != g.ref || l.Alt != g.alt {
			t.Errorf("row %d differs: varquery %s, gtcount %s:%d:%s:%s",
				i, l, g.chrom, g.pos, g.ref, g.alt)
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

// TestSampleFileForms pins that --sample takes a file as well as a name, the way
// --variant does.
func TestSampleFileForms(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")

	// A whitespace/newline list with comments.
	list := writeFile(t, "subjects.txt", "# subjects\nS1\nS2\n")
	got := parseNames(t, list)
	if strings.Join(got, ",") != "S1,S2" {
		t.Errorf("sample list read as %v", got)
	}

	// A VCF: its header roster, not its data lines. Read as a name list, a VCF would
	// silently yield thousands of bogus samples.
	got = parseNames(t, "testdata/coverage.vcf")
	if strings.Join(got, ",") != "S1,S2" {
		t.Errorf("VCF header roster read as %v", got)
	}

	// Names dedupe, first occurrence winning, so a repeat cannot become two columns.
	set, err := parseSampleArgs(context.Background(), []string{"S2", list, "S2"})
	if err != nil {
		t.Fatal(err)
	}
	if strings.Join(set.Names, ",") != "S2,S1" {
		t.Errorf("want S2,S1 (first-seen order, deduped), got %v", set.Names)
	}

	// End to end, the file selects the same rows the name does.
	viaFile := dataRowsOnly(runVcf(t, "vcf-varquery", "--sample",
		writeFile(t, "one.txt", "S1\n"), "--variant", "chr1:500:A:T", base))
	viaName := dataRowsOnly(runVcf(t, "vcf-varquery", "--sample", "S1",
		"--variant", "chr1:500:A:T", base))
	if viaFile != viaName {
		t.Errorf("a sample file and a sample name disagree:\n%s\n%s", viaFile, viaName)
	}
}

// TestUnknownSampleIsAnError pins that a sample the source lacks fails loudly.
//
// A store answers an unknown sample with silence, so without this a typo -- in a
// name, or in the path of a sample file, which is then taken as a name -- looks
// exactly like a subject that genuinely carries nothing.
func TestUnknownSampleIsAnError(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	for _, bad := range []string{"S99", "/nonexistent/subjects.txt"} {
		err := runVcfErr(t, "vcf-varquery", "--sample", bad, "--variant", "chr1:100:A:G", base)
		if err == nil {
			t.Errorf("%q should be rejected", bad)
			continue
		}
		if !strings.Contains(err.Error(), "not in this source") {
			t.Errorf("%q: unhelpful error %v", bad, err)
		}
	}
}

// TestDosageColumn pins --dosage, which is what PGS and GReX tooling consumes.
func TestDosageColumn(t *testing.T) {
	base := convert(t, "testdata/multiallelic.vcf")
	out := dataRowsOnly(runVcf(t, "vcf-varquery", "--variant", "chr1:200:C:T",
		"--variant", "chr1:100:A:G", "--hom-ref", "--dosage", base))

	header := strings.Split(out, "\n")[0]
	if !strings.HasSuffix(header, "\tdosage") {
		t.Fatalf("dosage should be appended, leaving the base layout alone: %q", header)
	}
	got := map[string]string{}
	for _, r := range strings.Split(out, "\n") {
		f := strings.Split(r, "\t")
		if len(f) != numCols+1 || f[colChrom] == "chrom" {
			continue
		}
		got[f[colSample]+"@"+f[colPos]+":"+f[colGT]] = f[numCols]
	}
	for key, want := range map[string]string{
		"S1@200:1/.": "1", // a split multiallelic: one copy of THIS alt
		"S2@200:0/1": "1",
		"S3@100:1/1": "2", // homozygous alt
		"S2@100:0/0": "0", // observed reference
	} {
		if got[key] != want {
			t.Errorf("dosage for %s = %q, want %q (all: %v)", key, got[key], want, got)
		}
	}
}

// parseNames is parseSampleArgs reduced to the names it resolved.
func parseNames(t *testing.T, arg string) []string {
	t.Helper()
	set, err := parseSampleArgs(context.Background(), []string{arg})
	if err != nil {
		t.Fatalf("parseSampleArgs(%s): %v", arg, err)
	}
	return set.Names
}
