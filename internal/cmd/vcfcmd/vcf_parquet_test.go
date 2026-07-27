package vcfcmd

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cgkit/internal/varstore"
)

// convert runs vcf-toparquet into a temp dir and returns the store base name.
func convert(t *testing.T, input string, extra ...string) string {
	t.Helper()
	base := filepath.Join(t.TempDir(), "store")
	args := append([]string{"vcf-toparquet", "--out", base}, extra...)
	runVcf(t, append(args, input)...)
	return base
}

// readCalls loads every call of a store back through the reading path, which is
// how output correctness is checked -- by re-parsing rather than diffing bytes.
func readCalls(t *testing.T, base string) []varstore.Call {
	t.Helper()
	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatalf("OpenParquet(%s): %v", base, err)
	}
	defer s.Close()

	samples, err := s.Samples()
	if err != nil {
		t.Fatalf("Samples: %v", err)
	}
	var calls []varstore.Call
	for _, name := range samples {
		c, err := s.Variants(name, nil, varstore.Gate{})
		if err != nil {
			t.Fatalf("Variants(%s): %v", name, err)
		}
		calls = append(calls, c...)
	}
	return calls
}

func TestVcfToParquetRequiresOut(t *testing.T) {
	err := runVcfErr(t, "vcf-toparquet", "testdata/multiallelic.vcf")
	if err == nil || !strings.Contains(err.Error(), "--out") {
		t.Errorf("expected an --out error, got %v", err)
	}
}

func TestVcfToParquetWritesAllThreeFiles(t *testing.T) {
	base := convert(t, "testdata/multiallelic.vcf")
	for _, p := range []string{
		varstore.CallsPath(base), varstore.SitesPath(base), varstore.RegionsPath(base),
	} {
		st, err := os.Stat(p)
		if err != nil {
			t.Errorf("missing store file %s: %v", p, err)
			continue
		}
		if st.Size() == 0 {
			t.Errorf("store file %s is empty", p)
		}
	}
}

// TestVcfToParquetSplitsMultiallelic pins the normalization contract: one row
// per ALT allele, the focal allele recoded to 1, other alternates masked to "."
// rather than to reference, and AD taken per allele rather than summed.
func TestVcfToParquetSplitsMultiallelic(t *testing.T) {
	base := convert(t, "testdata/multiallelic.vcf")
	calls := readCalls(t, base)

	got := map[string]varstore.Call{}
	for _, c := range calls {
		got[c.SampleID+"@"+c.Locus().String()] = c
	}

	// S1 is 1/2 at the multiallelic site: a carrier of BOTH alleles, appearing
	// once per split row, with the other allele masked out.
	t1, ok := got["S1@1:200:C:T"]
	if !ok {
		t.Fatalf("S1 missing from the C>T split row; got keys %v", keys(got))
	}
	if t1.GT != "1/." {
		t.Errorf("S1 C>T gt = %q, want %q", t1.GT, "1/.")
	}
	if t1.ADRef != 5 || t1.ADAlt != 12 {
		t.Errorf("S1 C>T ad = %d,%d, want 5,12", t1.ADRef, t1.ADAlt)
	}
	g1, ok := got["S1@1:200:C:G"]
	if !ok {
		t.Fatal("S1 missing from the C>G split row")
	}
	if g1.GT != "./1" {
		t.Errorf("S1 C>G gt = %q, want %q", g1.GT, "./1")
	}
	if g1.ADRef != 5 || g1.ADAlt != 9 {
		t.Errorf("S1 C>G ad = %d,%d, want 5,9 (per-allele, not summed)", g1.ADRef, g1.ADAlt)
	}

	// S4 is 2/2: a homozygous carrier of the second alt only.
	if c, ok := got["S4@1:200:C:G"]; !ok || c.GT != "1/1" {
		t.Errorf("S4 C>G = %+v, want gt 1/1", c)
	}
	if _, ok := got["S4@1:200:C:T"]; ok {
		t.Error("S4 should not appear as a carrier of C>T")
	}

	// A site where nobody carries anything produces no calls at all...
	if _, ok := got["S1@1:300:G:A"]; ok {
		t.Error("monomorphic site should produce no calls")
	}
}

// TestVcfToParquetSitesKeepsMonomorphic is the reason the sites file exists:
// a position everyone is reference at still has to be recorded as interrogated.
func TestVcfToParquetSitesKeepsMonomorphic(t *testing.T) {
	base := convert(t, "testdata/multiallelic.vcf")
	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	// chr1:300 has zero carriers, so it is absent from the calls entirely.
	// Classify must still call it interrogated, not not-assayed.
	states, err := s.Classify(varstore.Locus{Chrom: "1", Pos: 300, Ref: "G", Alt: "A"}, varstore.Gate{})
	if err != nil {
		t.Fatalf("Classify: %v", err)
	}
	for _, st := range states {
		if st.State != varstore.StateNonCarrier {
			t.Errorf("%s at monomorphic site = %s, want non_carrier", st.SampleID, st.State)
		}
	}
}

// TestVcfToParquetNoCallBreaksRun pins that depth alone does not make a site
// callable: S1 has DP 40 at chr1:400 but a "./." genotype, so its run must break.
func TestVcfToParquetNoCallBreaksRun(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	states, err := s.Classify(varstore.Locus{Chrom: "1", Pos: 400, Ref: "T", Alt: "C"}, varstore.Gate{})
	if err != nil {
		t.Fatalf("Classify: %v", err)
	}
	byName := map[string]varstore.State{}
	for _, st := range states {
		byName[st.SampleID] = st.State
	}
	if byName["S1"] != varstore.StateNotAssayed {
		t.Errorf("S1 at a no-call site = %s, want not_assayed (DP 40 but GT ./.)", byName["S1"])
	}
	if byName["S2"] != varstore.StateNonCarrier {
		t.Errorf("S2 at the same site = %s, want non_carrier", byName["S2"])
	}
}

// TestVcfToParquetLowDepthIsNotAssayed pins that a covered-looking sample below
// the threshold is not silently counted as a non-carrier.
func TestVcfToParquetLowDepthIsNotAssayed(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	states, err := s.Classify(varstore.Locus{Chrom: "1", Pos: 200, Ref: "C", Alt: "T"}, varstore.Gate{})
	if err != nil {
		t.Fatalf("Classify: %v", err)
	}
	for _, st := range states {
		if st.SampleID == "S2" && st.State != varstore.StateNotAssayed {
			t.Errorf("S2 at DP 3 (--min-dp 10) = %s, want not_assayed", st.State)
		}
	}
}

// TestVcfToParquetNoCallableRefuses pins that a source without DP is rejected
// rather than quietly producing a store that cannot classify.
func TestVcfToParquetNoCallableRefuses(t *testing.T) {
	err := runVcfErr(t, "vcf-toparquet", "--out", filepath.Join(t.TempDir(), "s"), "testdata/sample.vcf")
	if err == nil || !strings.Contains(err.Error(), "--no-callable") {
		t.Errorf("expected a --no-callable error for a DP-less VCF, got %v", err)
	}
}

// TestVarQueryClassifyRefusesWithoutRegions pins the loud-failure contract:
// an incomplete store must error, never report everyone as a non-carrier.
func TestVarQueryClassifyRefusesWithoutRegions(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	if err := os.Remove(varstore.RegionsPath(base)); err != nil {
		t.Fatal(err)
	}
	err := runVcfErr(t, "vcf-varquery", "--variant", "1:100:A:G", "--classify", base)
	if err == nil {
		t.Fatal("expected an error when the regions file is missing")
	}
	if !strings.Contains(err.Error(), "non-carrier from not-assayed") {
		t.Errorf("error should name the classification limit, got %v", err)
	}
}

// TestVarQueryBackendsAgree is the core equivalence claim: the same question
// against a VCF and against a store converted from it must give one answer.
func TestVarQueryBackendsAgree(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	for _, v := range []string{"1:100:A:G", "1:200:C:T", "1:400:T:C", "1:500:A:T"} {
		fromVcf := runVcf(t, "vcf-varquery", "--variant", v, "--classify", "--min-dp", "10", "testdata/coverage.vcf")
		fromPq := runVcf(t, "vcf-varquery", "--variant", v, "--classify", "--min-dp", "10", base)
		if dataRowsOnly(fromVcf) != dataRowsOnly(fromPq) {
			t.Errorf("backends disagree at %s\n vcf:\n%s\n parquet:\n%s", v,
				dataRowsOnly(fromVcf), dataRowsOnly(fromPq))
		}
	}
}

func TestVarQueryRejectsBadFormatForMode(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	if err := runVcfErr(t, "vcf-varquery", "--variant", "1:100:A:G", "--format", "vcf", base); err == nil {
		t.Error("--format vcf should be rejected in --variant mode")
	}
	if err := runVcfErr(t, "vcf-varquery", "--sample", "S1", "--format", "list", base); err == nil {
		t.Error("--format list should be rejected in --sample mode")
	}
	if err := runVcfErr(t, "vcf-varquery", "--sample", "S1", "--variant", "1:100:A:G", base); err == nil {
		t.Error("--sample and --variant together should be rejected")
	}
}

// dataRowsOnly drops the ## provenance lines, which carry a timestamp.
func dataRowsOnly(s string) string {
	var out []string
	for _, l := range strings.Split(s, "\n") {
		if !strings.HasPrefix(l, "##") && l != "" {
			out = append(out, l)
		}
	}
	return strings.Join(out, "\n")
}

func keys(m map[string]varstore.Call) []string {
	out := make([]string, 0, len(m))
	for k := range m {
		out = append(out, k)
	}
	return out
}
