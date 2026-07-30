package vcfcmd

import (
	"os"
	"strings"
	"testing"

	"github.com/compgenlab/cgkit/internal/varstore"
)

// Tests for --hom-ref, which adds reference calls to either query mode.
//
// The recurring concern is that a reference call is an assertion about data, not
// a default: it may only be reported where the source actually observed the
// sample to be all-reference. So most of these pin what must NOT appear.

// Column positions in the shared tabular layout, which both query modes use:
//
//	chrom pos ref alt sample gt dp min_dp ad_ref ad_alt gq
const (
	colChrom = iota
	colPos
	colRef
	colAlt
	colSample
	colGT
	colDP
	colMinDP
	colADRef
	colADAlt
	colGQ
	numCols
)

// identityRows keeps the columns that identify a row -- through gt, dropping the
// quality columns.
//
// Those are deliberately excluded: a Parquet store keeps only alternate
// genotypes, so a reference call recovered from one has no DP/AD/GQ to report
// where the same query against a VCF does, and min_dp is a threshold on one
// backend and an exact depth on the other. That asymmetry is real and documented;
// what must not differ is WHICH rows appear.
func identityRows(s string) string {
	var out []string
	for _, l := range strings.Split(dataRowsOnly(s), "\n") {
		f := strings.Split(l, "\t")
		if len(f) > colDP {
			f = f[:colDP]
		}
		out = append(out, strings.Join(f, "\t"))
	}
	return strings.Join(out, "\n")
}

// TestHomRefBackendsAgree extends the core equivalence claim to reference calls:
// a VCF and a store converted from it must name the same samples as 0/0.
func TestHomRefBackendsAgree(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	for _, v := range []string{"1:100:A:G", "1:200:C:T", "1:400:T:C", "1:500:A:T", "2:100:G:C"} {
		fromVcf := runVcf(t, "vcf-varquery", "--variant", v, "--hom-ref",
			"--min-dp", "10", "testdata/coverage.vcf")
		fromPq := runVcf(t, "vcf-varquery", "--variant", v, "--hom-ref", "--min-dp", "10", base)
		if identityRows(fromVcf) != identityRows(fromPq) {
			t.Errorf("backends disagree at %s\n vcf:\n%s\n parquet:\n%s", v,
				identityRows(fromVcf), identityRows(fromPq))
		}
	}
}

func TestHomRefSampleModeBackendsAgree(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	for _, s := range []string{"S1", "S2"} {
		fromVcf := runVcf(t, "vcf-varquery", "--sample", s, "--hom-ref",
			"--min-dp", "10", "testdata/coverage.vcf")
		fromPq := runVcf(t, "vcf-varquery", "--sample", s, "--hom-ref", "--min-dp", "10", base)
		if identityRows(fromVcf) != identityRows(fromPq) {
			t.Errorf("backends disagree for %s\n vcf:\n%s\n parquet:\n%s", s,
				identityRows(fromVcf), identityRows(fromPq))
		}
	}
}

// TestHomRefIsStricterThanNonCarrier is the central semantic test.
//
// At the multiallelic record chr1:200 C>T,G sample S2 is 0/1: not a carrier of
// G, so --classify calls it non_carrier there -- but its genotype is not
// reference either, and reporting 0/0 for it would be a genotype the source
// never contained. S3 is the only genuinely all-reference sample.
func TestHomRefIsStricterThanNonCarrier(t *testing.T) {
	base := convert(t, "testdata/multiallelic.vcf")
	for _, in := range []string{"testdata/multiallelic.vcf", base} {
		homRef := homRefSamples(t, runVcf(t, "vcf-varquery",
			"--variant", "chr1:200:C:G", "--hom-ref", in))
		if len(homRef) != 1 || homRef[0] != "S3" {
			t.Errorf("%s: hom-ref samples at chr1:200 C>G = %v, want [S3] only\n"+
				"S2 is 0/1 and S4 is 2/2 there: neither carries G, but neither is reference",
				in, homRef)
		}

		// The same query under --classify must still call S2 a non-carrier, so
		// this is a real distinction between the two flags rather than a change
		// to what non_carrier means.
		classify := runVcf(t, "vcf-varquery", "--variant", "chr1:200:C:G", "--classify", in)
		if !strings.Contains(dataRowsOnly(classify), "S2\tnon_carrier") {
			t.Errorf("%s: --classify should still report S2 as non_carrier at chr1:200 C>G:\n%s",
				in, dataRowsOnly(classify))
		}
	}
}

// TestHomRefSampleModeSkipsOtherAlternates is the same rule from the sample
// side: S2 carries T at chr1:200, so neither split site of that record may be
// reported as reference for it.
func TestHomRefSampleModeSkipsOtherAlternates(t *testing.T) {
	base := convert(t, "testdata/multiallelic.vcf")
	for _, in := range []string{"testdata/multiallelic.vcf", base} {
		got := runVcf(t, "vcf-varquery", "--sample", "S2", "--hom-ref", in)
		for _, row := range strings.Split(dataRowsOnly(got), "\n") {
			f := strings.Split(row, "\t")
			if len(f) < numCols || f[colSample] != "S2" {
				continue
			}
			if f[colPos] == "200" && f[colGT] == varstore.HomRefGT {
				t.Errorf("%s: S2 reported as 0/0 at the multiallelic record it carries T at:\n%s",
					in, row)
			}
		}
	}
}

// TestHomRefUnionKeepsCarriers pins that --hom-ref adds to the output rather
// than replacing it, and that the gt column tells the two kinds apart.
func TestHomRefUnionKeepsCarriers(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	got := dataRowsOnly(runVcf(t, "vcf-varquery", "--sample", "S1", "--hom-ref", "--min-dp", "10", base))
	if !strings.Contains(got, "\t500\tA\tT\tS1\t0/1\t") {
		t.Errorf("the carried variant at chr1:500 should still be reported:\n%s", got)
	}
	if !strings.Contains(got, "\t300\tG\tA\tS1\t0/0\t") {
		t.Errorf("the reference call at chr1:300 should be reported:\n%s", got)
	}

	// Without the flag, only the carried variant appears.
	plain := dataRowsOnly(runVcf(t, "vcf-varquery", "--sample", "S1", "--min-dp", "10", base))
	if strings.Contains(plain, varstore.HomRefGT) {
		t.Errorf("reference calls must not appear without --hom-ref:\n%s", plain)
	}
}

// TestHomRefNeedsAnObservation pins that a no-call and a below-depth call are
// both excluded. S1 is "./." at DP 40 at chr1:400 -- a declined call, not a
// reference one -- and S2 is 0/0 at DP 3 at chr1:200, which --min-dp 10 rejects.
func TestHomRefNeedsAnObservation(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	for _, in := range []string{"testdata/coverage.vcf", base} {
		at400 := homRefSamples(t, runVcf(t, "vcf-varquery",
			"--variant", "1:400:T:C", "--hom-ref", "--min-dp", "10", in))
		if contains(at400, "S1") {
			t.Errorf("%s: S1 is ./. at chr1:400 (DP 40); a declined call is not a "+
				"reference call, got %v", in, at400)
		}
		at200 := homRefSamples(t, runVcf(t, "vcf-varquery",
			"--variant", "1:200:C:T", "--hom-ref", "--min-dp", "10", in))
		if contains(at200, "S2") {
			t.Errorf("%s: S2 is 0/0 at DP 3 at chr1:200, below --min-dp 10, got %v", in, at200)
		}
	}
}

// TestHomRefGateAppliesToVcfBackend pins that the gate acts on the reference
// call itself, so the same VCF gives different answers at different thresholds.
// (A store cannot vary this after the fact: its runs baked in the conversion
// --min-dp, which is what the verbose mismatch warning is about.)
func TestHomRefGateAppliesToVcfBackend(t *testing.T) {
	ungated := homRefSamples(t, runVcf(t, "vcf-varquery",
		"--variant", "1:200:C:T", "--hom-ref", "testdata/coverage.vcf"))
	if !contains(ungated, "S2") {
		t.Errorf("ungated, S2's 0/0 at DP 3 should be reported, got %v", ungated)
	}
	gated := homRefSamples(t, runVcf(t, "vcf-varquery",
		"--variant", "1:200:C:T", "--hom-ref", "--min-gq", "50", "testdata/coverage.vcf"))
	if contains(gated, "S2") {
		t.Errorf("--min-gq 50 should exclude S2's GQ 20 reference call, got %v", gated)
	}
}

// TestHomRefOffCatalogReportsNothing pins that the sites catalog still bounds
// what can be answered. chr1:250 lies between records and is bracketed by both
// samples' runs, so a store reading runs as coverage would invent two reference
// calls there.
func TestHomRefOffCatalogReportsNothing(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	for _, in := range []string{"testdata/coverage.vcf", base} {
		got := homRefSamples(t, runVcf(t, "vcf-varquery",
			"--variant", "1:250:A:G", "--hom-ref", in))
		if len(got) != 0 {
			t.Errorf("%s: an off-catalog locus must yield no reference calls, got %v", in, got)
		}
	}
}

// TestHomRefRefusesIncompleteStore pins the loud-failure contract: a store that
// cannot tell a reference call from an unassayed position must error in both
// modes rather than silently report nobody as reference.
func TestHomRefRefusesIncompleteStore(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	if err := os.Remove(varstore.RegionsPath(base)); err != nil {
		t.Fatal(err)
	}
	for _, args := range [][]string{
		{"vcf-varquery", "--variant", "1:100:A:G", "--hom-ref", base},
		{"vcf-varquery", "--sample", "S1", "--hom-ref", base},
	} {
		err := runVcfErr(t, args...)
		if err == nil {
			t.Fatalf("%v: expected an error when the regions file is missing", args)
		}
		if !strings.Contains(err.Error(), "non-carrier from not-assayed") {
			t.Errorf("%v: error should name the classification limit, got %v", args, err)
		}
	}
}

// TestHomRefRejectedWhereItCannotBeRead pins the two combinations that would
// produce output nobody can interpret.
func TestHomRefRejectedWhereItCannotBeRead(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	if err := runVcfErr(t, "vcf-varquery", "--variant", "1:100:A:G",
		"--hom-ref", "--classify", base); err == nil {
		t.Error("--hom-ref with --classify should be rejected; --classify already resolves every sample")
	}
	if err := runVcfErr(t, "vcf-varquery", "--variant", "1:100:A:G",
		"--hom-ref", "--format", "list", base); err == nil {
		t.Error("--hom-ref with --format list should be rejected; bare ids cannot say which are carriers")
	}
}

// TestHomRefStoreRowsCarryNoQuality pins the documented asymmetry, so it stays a
// known limit rather than becoming a surprise: a reference call from a store has
// no DP/AD/GQ, while a carrier row alongside it still does.
func TestHomRefStoreRowsCarryNoQuality(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	for _, row := range strings.Split(dataRowsOnly(runVcf(t, "vcf-varquery",
		"--sample", "S1", "--hom-ref", "--min-dp", "10", base)), "\n") {
		f := strings.Split(row, "\t")
		if len(f) != numCols || f[colSample] != "S1" {
			continue
		}
		gt, dp, gq := f[colGT], f[colDP], f[colGQ]
		if gt == varstore.HomRefGT && (dp != "." || gq != ".") {
			t.Errorf("a store cannot know a reference call's quality, got dp=%s gq=%s:\n%s",
				dp, gq, row)
		}
		if gt == "0/1" && dp == "." {
			t.Errorf("a carrier row must keep its depth:\n%s", row)
		}
	}
}

// TestBothModesShareOneLayout pins the output contract: --sample and --variant
// emit identical columns, so their output can be concatenated, cut and sorted the
// same way. They used to differ -- --sample led with the sample column.
func TestBothModesShareOneLayout(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	want := strings.Join([]string{
		"chrom", "pos", "ref", "alt", "sample", "gt", "dp", "min_dp", "ad_ref", "ad_alt", "gq",
	}, "\t")
	for _, args := range [][]string{
		{"vcf-varquery", "--sample", "S1", base},
		{"vcf-varquery", "--variant", "chr1:100:A:G", base},
		{"vcf-varquery", "--sample", "S1", "--hom-ref", base},
		{"vcf-varquery", "--variant", "chr1:100:A:G", "--hom-ref", base},
	} {
		header := strings.Split(dataRowsOnly(runVcf(t, args...)), "\n")[0]
		if header != want {
			t.Errorf("%v\n header = %q\n   want = %q", args[1:], header, want)
		}
	}
}

// TestMinDPReportsTheVouchedFloor pins what min_dp means: the tightest lower
// bound on depth the backend can vouch for.
//
// This is the whole reason the column exists. A reference call recovered from a
// store has no recorded depth, but it came from a run built at the conversion
// --min-dp -- so the store can still say the site was covered well enough to
// call, which a bare "." would throw away.
func TestMinDPReportsTheVouchedFloor(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf") // --min-dp defaults to 10

	rows := map[string][]string{}
	for _, l := range strings.Split(dataRowsOnly(runVcf(t, "vcf-varquery",
		"--sample", "S1", "--hom-ref", "--min-dp", "10", base)), "\n") {
		if f := strings.Split(l, "\t"); len(f) == numCols && f[colChrom] != "chrom" {
			rows[f[colChrom]+":"+f[colPos]] = f
		}
	}

	// A reconstructed reference call: no depth of its own, but the run vouches
	// for the conversion threshold.
	ref, ok := rows["chr1:300"]
	if !ok {
		t.Fatalf("expected a reference call at chr1:300; got %v", rows)
	}
	if ref[colDP] != "." || ref[colMinDP] != "10" {
		t.Errorf("store reference call: dp=%s min_dp=%s, want dp=. min_dp=10 "+
			"(no recorded depth, but callable at the conversion threshold)",
			ref[colDP], ref[colMinDP])
	}

	// A carrier records its own depth, which is its own best bound.
	carrier, ok := rows["chr1:500"]
	if !ok {
		t.Fatalf("expected a carrier at chr1:500; got %v", rows)
	}
	if carrier[colDP] != "30" || carrier[colMinDP] != "30" {
		t.Errorf("carrier: dp=%s min_dp=%s, want both 30", carrier[colDP], carrier[colMinDP])
	}

	// A VCF backend knows every depth exactly, so min_dp is never a threshold.
	for _, l := range strings.Split(dataRowsOnly(runVcf(t, "vcf-varquery",
		"--sample", "S1", "--hom-ref", "--min-dp", "10", "testdata/coverage.vcf")), "\n") {
		f := strings.Split(l, "\t")
		if len(f) != numCols || f[colChrom] == "chrom" {
			continue
		}
		if f[colDP] != f[colMinDP] {
			t.Errorf("vcf backend: dp=%s but min_dp=%s; an exact depth is its own bound\n%s",
				f[colDP], f[colMinDP], l)
		}
	}
}

// TestHomRefRegionBounds pins that --region restricts the catalog walk, which is
// the only thing keeping a whole-genome --sample query bounded.
func TestHomRefRegionBounds(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	got := dataRowsOnly(runVcf(t, "vcf-varquery", "--sample", "S2", "--hom-ref",
		"--min-dp", "10", "--region", "chr1:250-450", base))
	for _, want := range []string{"\t300\t", "\t400\t"} {
		if !strings.Contains(got, want) {
			t.Errorf("expected a row at %s within the region:\n%s", want, got)
		}
	}
	for _, unwanted := range []string{"\t100\t", "\t200\t", "\t500\t"} {
		if strings.Contains(got, unwanted) {
			t.Errorf("row at %s is outside --region chr1:250-450:\n%s", unwanted, got)
		}
	}
}

// TestVarQueryVcfOutputIsSorted pins that --format vcf emits coordinate-sorted
// records. Loci are gathered per sample, so without an explicit sort the second
// sample's private loci all land after the first sample's -- and an unsorted VCF
// cannot be indexed. --hom-ref makes that the normal case rather than a rare one.
func TestVarQueryVcfOutputIsSorted(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	out := runVcf(t, "vcf-varquery", "--sample", "S1", "--sample", "S2",
		"--hom-ref", "--min-dp", "10", "--format", "vcf", base)

	var last string
	var lastPos int64
	for _, l := range strings.Split(out, "\n") {
		if l == "" || strings.HasPrefix(l, "#") {
			continue
		}
		f := strings.Split(l, "\t")
		pos := atoiTest(t, f[1])
		if f[0] == last && pos < lastPos {
			t.Errorf("VCF output is not coordinate sorted at %s:%d (previous %d)\n%s",
				f[0], pos, lastPos, out)
		}
		if f[0] != last {
			if last != "" && f[0] < last {
				t.Errorf("chromosomes are out of order: %s after %s\n%s", f[0], last, out)
			}
			last = f[0]
		}
		lastPos = pos
	}
}

// homRefSamples extracts the sample ids reported as reference calls from a
// --variant query, whose sample column is fifth and gt column sixth.
func homRefSamples(t *testing.T, out string) []string {
	t.Helper()
	var got []string
	for _, l := range strings.Split(dataRowsOnly(out), "\n") {
		f := strings.Split(l, "\t")
		if len(f) < numCols || f[colChrom] == "chrom" {
			continue
		}
		if f[colGT] == varstore.HomRefGT {
			got = append(got, f[colSample])
		}
	}
	return got
}

func contains(haystack []string, want string) bool {
	for _, h := range haystack {
		if h == want {
			return true
		}
	}
	return false
}

func atoiTest(t *testing.T, s string) int64 {
	t.Helper()
	var n int64
	for _, c := range s {
		if c < '0' || c > '9' {
			t.Fatalf("not a position: %q", s)
		}
		n = n*10 + int64(c-'0')
	}
	return n
}
