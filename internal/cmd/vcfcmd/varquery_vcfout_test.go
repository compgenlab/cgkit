package vcfcmd

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// Tests for --format vcf: the genotype matrix, its recomputed INFO fields, and
// the record ordering an indexable VCF requires.

// vcfRecords returns the data lines of a VCF, dropping ## and #CHROM.
func vcfRecords(out string) []string {
	var rows []string
	for _, l := range strings.Split(out, "\n") {
		if l == "" || strings.HasPrefix(l, "#") {
			continue
		}
		rows = append(rows, l)
	}
	return rows
}

// info pulls one INFO key out of a VCF record.
func info(t *testing.T, record, key string) string {
	t.Helper()
	f := strings.Split(record, "\t")
	if len(f) < 8 {
		t.Fatalf("record has %d columns: %q", len(f), record)
	}
	for _, kv := range strings.Split(f[7], ";") {
		if k, v, ok := strings.Cut(kv, "="); ok && k == key {
			return v
		}
	}
	t.Fatalf("no %s in INFO %q", key, f[7])
	return ""
}

// TestVcfOutRecordsAreInContigOrder is the regression test for a bug that shipped.
//
// The old writer sorted loci with a string comparison on chrom, so chr10 came out
// before chr2. Harmless on a two-contig fixture and wrong for any real store,
// because an indexable VCF needs records grouped in the header's contig order.
// Streaming in the store's own write order gives that with no comparator at all.
func TestVcfOutRecordsAreInContigOrder(t *testing.T) {
	base := convert(t, "testdata/contigs.vcf")
	for _, in := range []string{"testdata/contigs.vcf", base} {
		rows := vcfRecords(runVcf(t, "vcf-varquery", "--variant", "chr2", "--variant", "chr10",
			"--format", "vcf", "--min-dp", "10", in))
		if len(rows) != 4 {
			t.Fatalf("%s: want 4 records, got %d: %v", in, len(rows), rows)
		}
		var order []string
		for _, r := range rows {
			order = append(order, strings.SplitN(r, "\t", 2)[0])
		}
		got := strings.Join(order, ",")
		if got != "chr2,chr2,chr10,chr10" {
			t.Errorf("%s: records in %s; want chr2,chr2,chr10,chr10 (contig order as written, "+
				"not lexicographic -- a string sort puts chr10 first and tabix cannot index that)",
				in, got)
		}
	}
}

// TestVcfOutImpliesHomRef pins that a genotype matrix reports observed reference
// calls as 0/0 rather than ./.. Without that every non-carrier is indistinguishable
// from an unassayed sample, and the file asserts far less than the data supports.
func TestVcfOutImpliesHomRef(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	rows := vcfRecords(runVcf(t, "vcf-varquery", "--variant", "chr1:100:A:G",
		"--format", "vcf", "--min-dp", "10", base))
	if len(rows) != 1 {
		t.Fatalf("want one record, got %v", rows)
	}
	// S1 is 0/0 at chr1:100 and S2 is 0/1. No --hom-ref was passed.
	f := strings.Split(rows[0], "\t")
	s1, s2 := f[len(f)-2], f[len(f)-1]
	if !strings.HasPrefix(s1, "0/0") {
		t.Errorf("S1 should be 0/0, not %q -- --format vcf must imply reference calls", s1)
	}
	if !strings.HasPrefix(s2, "0/1") {
		t.Errorf("S2 should be 0/1, got %q", s2)
	}
}

// TestVcfOutInfoMatchesStoreWhenComplete is the consistency check that falls out
// of recomputing AC/AN: with every sample included and nothing below the depth
// threshold, the recomputed values must equal the store's own catalog.
func TestVcfOutInfoMatchesStoreWhenComplete(t *testing.T) {
	base := convert(t, "testdata/contigs.vcf")
	rows := vcfRecords(runVcf(t, "vcf-varquery", "--variant", "chr2",
		"--format", "vcf", "--min-dp", "10", base))
	if len(rows) != 2 {
		t.Fatalf("want 2 records, got %v", rows)
	}
	// chr2:200 -- S1 is 0/0, S2 is 1/1. AC counts alleles, so 2; AN is 4.
	rec := rows[1]
	for key, want := range map[string]string{"AC": "2", "AN": "4", "AF": "0.5", "NS": "2", "nhomalt": "1"} {
		if got := info(t, rec, key); got != want {
			t.Errorf("chr2:200 %s = %s, want %s (record: %s)", key, got, want, rec)
		}
	}
}

// TestVcfOutInfoIsRecomputedOverTheSubset is the reason the values are recomputed
// rather than copied from sites.parquet: over a sample subset the store's counts
// are for the wrong denominator.
func TestVcfOutInfoIsRecomputedOverTheSubset(t *testing.T) {
	base := convert(t, "testdata/contigs.vcf")

	all := vcfRecords(runVcf(t, "vcf-varquery", "--variant", "chr2",
		"--format", "vcf", "--min-dp", "10", base))
	// S2 alone at chr2:200, where it is 1/1: AC 2 of AN 2, so AF is 1.0 -- not the
	// 0.5 the whole-cohort catalog reports.
	sub := vcfRecords(runVcf(t, "vcf-varquery", "--variant", "chr2", "--sample", "S2",
		"--format", "vcf", "--min-dp", "10", base))
	if len(all) != 2 || len(sub) != 2 {
		t.Fatalf("want 2 records each, got %d and %d", len(all), len(sub))
	}
	if got := info(t, sub[1], "AN"); got != "2" {
		t.Errorf("one-sample AN = %s, want 2", got)
	}
	if got := info(t, sub[1], "AF"); got != "1" {
		t.Errorf("one-sample AF = %s, want 1 (S2 is 1/1); the cohort value is %s",
			got, info(t, all[1], "AF"))
	}
	// And the column count follows the subset.
	if n := len(strings.Split(sub[1], "\t")); n != 10 {
		t.Errorf("a one-sample matrix should have 10 columns, got %d: %s", n, sub[1])
	}
}

// TestVcfOutBackendsAgree pins that the matrix is the same from a VCF and from a
// store converted from it, INFO included.
func TestVcfOutBackendsAgree(t *testing.T) {
	base := convert(t, "testdata/contigs.vcf")
	a := vcfRecords(runVcf(t, "vcf-varquery", "--variant", "chr2", "--variant", "chr10",
		"--format", "vcf", "--min-dp", "10", "testdata/contigs.vcf"))
	b := vcfRecords(runVcf(t, "vcf-varquery", "--variant", "chr2", "--variant", "chr10",
		"--format", "vcf", "--min-dp", "10", base))
	// The GT columns carry real DP/GQ from a VCF and none from a store, so compare
	// the locus and INFO columns -- which is where the recomputation lives.
	strip := func(rows []string) []string {
		out := make([]string, 0, len(rows))
		for _, r := range rows {
			f := strings.Split(r, "\t")
			out = append(out, strings.Join(f[:8], "\t"))
		}
		return out
	}
	if strings.Join(strip(a), "|") != strings.Join(strip(b), "|") {
		t.Errorf("backends disagree\n vcf:\n  %s\n parquet:\n  %s",
			strings.Join(strip(a), "\n  "), strings.Join(strip(b), "\n  "))
	}
}

// TestVcfOutBgzipAndIndex pins that --format vcf goes through cghts's VcfWriter
// rather than hand-written text, which is what makes bgzip and --tbi available.
//
// It used to write plain text with fmt.Fprintln, so a genotype matrix over a real
// cohort was uncompressed and unindexable -- close to useless as PGS input.
func TestVcfOutBgzipAndIndex(t *testing.T) {
	base := convert(t, "testdata/contigs.vcf")
	path := filepath.Join(t.TempDir(), "matrix.vcf.gz")
	runVcf(t, "vcf-varquery", "--variant", "chr2", "--variant", "chr10",
		"--min-dp", "10", "--format", "vcf", "-o", path, "--tbi", base)

	// BGZF sets the gzip FEXTRA flag; stdlib gzip does not, so this distinguishes a
	// regression to plain gzip.
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if len(raw) < 4 || raw[3]&0x04 == 0 {
		t.Errorf("output is gzip but not BGZF (flags=0x%02x); tabix cannot index it", raw[3])
	}
	if _, err := os.Stat(path + ".tbi"); err != nil {
		t.Errorf("--tbi did not write an index: %v", err)
	}

	// The written file must be readable as a store in its own right, which exercises
	// the header, the records and the genotypes together.
	rows := tsvDataRows(dataRowsOnly(runVcf(t, "vcf-varquery",
		"--variant", "chr10:200:T:C", "--hom-ref", "--min-dp", "10", path)))
	if len(rows) != 2 {
		t.Fatalf("round-trip gave %d rows, want 2: %v", len(rows), rows)
	}
	if !strings.Contains(rows[0], "\tS1\t0/0\t") || !strings.Contains(rows[1], "\tS2\t0/1\t") {
		t.Errorf("round-tripped genotypes wrong: %v", rows)
	}
}

// TestVcfOutTbiNeedsBgzipName pins that --tbi refuses a name it cannot index
// rather than writing an index nothing can use.
func TestVcfOutTbiNeedsBgzipName(t *testing.T) {
	base := convert(t, "testdata/contigs.vcf")
	err := runVcfErr(t, "vcf-varquery", "--variant", "chr2", "--format", "vcf",
		"-o", filepath.Join(t.TempDir(), "plain.vcf"), "--tbi", base)
	if err == nil {
		t.Fatal("--tbi with a non-bgzip output name should be refused")
	}
	if !strings.Contains(err.Error(), ".gz") {
		t.Errorf("the error should name the requirement, got %v", err)
	}
}
