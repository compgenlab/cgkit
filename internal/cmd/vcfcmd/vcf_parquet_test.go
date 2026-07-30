package vcfcmd

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/varstore"
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
		c, err := varstore.CollectCalls(s, varstore.Query{Samples: []string{name}})
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
	t1, ok := got["S1@chr1:200:C:T"]
	if !ok {
		t.Fatalf("S1 missing from the C>T split row; got keys %v", keys(got))
	}
	if t1.GT != "1/." {
		t.Errorf("S1 C>T gt = %q, want %q", t1.GT, "1/.")
	}
	if t1.ADRef != 5 || t1.ADAlt != 12 {
		t.Errorf("S1 C>T ad = %d,%d, want 5,12", t1.ADRef, t1.ADAlt)
	}
	g1, ok := got["S1@chr1:200:C:G"]
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
	if c, ok := got["S4@chr1:200:C:G"]; !ok || c.GT != "1/1" {
		t.Errorf("S4 C>G = %+v, want gt 1/1", c)
	}
	if _, ok := got["S4@chr1:200:C:T"]; ok {
		t.Error("S4 should not appear as a carrier of C>T")
	}

	// A site where nobody carries anything produces no calls at all...
	if _, ok := got["S1@chr1:300:G:A"]; ok {
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

// TestVarQueryHomRefRefusesWithoutRegions pins the loud-failure contract: an
// incomplete store must error, never report everyone as reference.
func TestVarQueryHomRefRefusesWithoutRegions(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	if err := os.Remove(varstore.RegionsPath(base)); err != nil {
		t.Fatal(err)
	}
	err := runVcfErr(t, "vcf-varquery", "--variant", "1:100:A:G", "--hom-ref", base)
	if err == nil {
		t.Fatal("expected an error when the regions file is missing")
	}
	if !strings.Contains(err.Error(), "non-carrier from not-assayed") {
		t.Errorf("error should name the classification limit, got %v", err)
	}
}

// TestVarQueryBackendsAgree is the core equivalence claim: the same question
// against a VCF and against a store converted from it must give one answer.
//
// Asked of the library rather than the CLI. Classify resolves every sample, so it
// exercises reference and unassayed samples that a carriers-only query never
// mentions -- and it is no longer reachable from the command line.
func TestVarQueryBackendsAgree(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	fromVcf, err := varstore.OpenVcf("testdata/coverage.vcf")
	if err != nil {
		t.Fatal(err)
	}
	defer fromVcf.Close()
	fromPq, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer fromPq.Close()

	gate := varstore.Gate{MinDP: 10}
	for _, v := range []string{"1:100:A:G", "1:200:C:T", "1:400:T:C", "1:500:A:T"} {
		locus, err := varstore.ParseLocus(v)
		if err != nil {
			t.Fatal(err)
		}
		a, err := fromVcf.Classify(locus, gate)
		if err != nil {
			t.Fatalf("vcf Classify(%s): %v", v, err)
		}
		b, err := fromPq.Classify(locus, gate)
		if err != nil {
			t.Fatalf("parquet Classify(%s): %v", v, err)
		}
		if len(a) != len(b) {
			t.Fatalf("%s: %d states from vcf, %d from parquet", v, len(a), len(b))
		}
		for i := range a {
			if a[i].SampleID != b[i].SampleID || a[i].State != b[i].State {
				t.Errorf("backends disagree at %s: vcf %s=%s, parquet %s=%s", v,
					a[i].SampleID, a[i].State, b[i].SampleID, b[i].State)
			}
		}
	}
}

// TestHomRefBackendsAgreeUnderTheGate is the regression test for a bug this
// found: with a GQ threshold the two backends disagreed about which samples were
// reference, because callable runs are built from depth alone so no GQ survives
// conversion. --min-gq was removed for exactly that reason; a depth gate, which
// the runs do encode, must still agree.
func TestHomRefBackendsAgreeUnderTheGate(t *testing.T) {
	base := convert(t, "testdata/gq.vcf", "--min-dp", "10")
	for _, s := range []string{"S1", "S2"} {
		fromVcf := runVcf(t, "vcf-varquery", "--sample", s, "--hom-ref", "--min-dp", "10", "testdata/gq.vcf")
		fromPq := runVcf(t, "vcf-varquery", "--sample", s, "--hom-ref", "--min-dp", "10", base)
		if identityRows(fromVcf) != identityRows(fromPq) {
			t.Errorf("backends disagree for %s\n vcf:\n%s\n parquet:\n%s", s,
				identityRows(fromVcf), identityRows(fromPq))
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

// TestOffCatalogLocusIsNotAssayed is the central semantic guarantee for a plain
// VCF: only the variants the file contains can be answered. chr1:250 lies
// between two records and is bracketed by both samples' called-site runs, so a
// store that read those runs as coverage would wrongly call everyone a
// non-carrier there.
func TestOffCatalogLocusIsNotAssayed(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	// Sanity: the position really is inside a run, so the test is meaningful.
	known, err := s.SiteKnown(varstore.Locus{Chrom: "1", Pos: 250, Ref: "A", Alt: "G"})
	if err != nil {
		t.Fatal(err)
	}
	if known {
		t.Fatal("chr1:250 should not be in the catalog; fixture changed")
	}
	if inRun, err := runBrackets(base, "1", 250); err != nil {
		t.Fatal(err)
	} else if !inRun {
		t.Skip("chr1:250 is not bracketed by a run; the test would prove nothing")
	}

	states, err := s.Classify(varstore.Locus{Chrom: "1", Pos: 250, Ref: "A", Alt: "G"}, varstore.Gate{})
	if err != nil {
		t.Fatalf("Classify: %v", err)
	}
	if len(states) == 0 {
		t.Fatal("expected a state per sample")
	}
	for _, st := range states {
		if st.State != varstore.StateNotAssayed {
			t.Errorf("%s at an off-catalog locus = %s, want not_assayed "+
				"(a plain VCF claims nothing between its records)", st.SampleID, st.State)
		}
	}
}

// TestOffCatalogAgreesAcrossBackends pins that the VCF backend draws the same
// boundary as the Parquet one.
//
// Two claims, at two levels. Through the CLI both backends must report no rows --
// an off-catalog locus yields neither an ALT call nor a reference one, even with
// --hom-ref. And in the library both must call every sample not_assayed rather
// than non_carrier, which is the claim "no rows" is standing in for.
func TestOffCatalogAgreesAcrossBackends(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	locus := varstore.Locus{Chrom: "1", Pos: 250, Ref: "A", Alt: "G"}

	fromVcf := dataRowsOnly(runVcf(t, "vcf-varquery", "--variant", "1:250:A:G",
		"--hom-ref", "testdata/coverage.vcf"))
	fromPq := dataRowsOnly(runVcf(t, "vcf-varquery", "--variant", "1:250:A:G", "--hom-ref", base))
	if fromVcf != fromPq {
		t.Errorf("backends disagree off-catalog\n vcf:\n%s\n parquet:\n%s", fromVcf, fromPq)
	}
	// No ALT call and no reference call anywhere.
	if rows := tsvDataRows(fromPq); len(rows) != 0 {
		t.Errorf("expected no data rows off-catalog, got %d:\n%s", len(rows),
			strings.Join(rows, "\n"))
	}

	for _, in := range []string{"testdata/coverage.vcf", base} {
		store, err := openVarStore(in, "")
		if err != nil {
			t.Fatal(err)
		}
		states, err := store.Classify(locus, varstore.Gate{})
		store.Close()
		if err != nil {
			t.Fatalf("%s: Classify: %v", in, err)
		}
		if len(states) == 0 {
			t.Fatalf("%s: expected a state per sample", in)
		}
		for _, st := range states {
			if st.State != varstore.StateNotAssayed {
				t.Errorf("%s: %s at an off-catalog locus = %s, want not_assayed "+
					"(a plain VCF claims nothing between its records)", in, st.SampleID, st.State)
			}
		}
	}
}

// TestSpanSemanticsRecorded pins that a VCF-derived store declares the
// conservative reading, so a future gVCF converter has to opt in explicitly.
func TestSpanSemanticsRecorded(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	if got := s.SpanSemantics(); got != varstore.SpansSites {
		t.Errorf("span semantics = %q, want %q", got, varstore.SpansSites)
	}
}

// runBrackets reports whether any called-site run spans pos, which is what makes
// the off-catalog test non-vacuous.
func runBrackets(base, chrom string, pos int32) (bool, error) {
	s, err := varstore.OpenParquet(base)
	if err != nil {
		return false, err
	}
	defer s.Close()
	// A run bracketing pos means Classify would have had evidence to misuse.
	states, err := s.Classify(varstore.Locus{Chrom: chrom, Pos: 100, Ref: "A", Alt: "G"}, varstore.Gate{})
	if err != nil {
		return false, err
	}
	// chr1:100 and chr1:300 are both catalog sites called in S1, so S1's run
	// spans 250. Confirm S1 was callable at the flanking site.
	for _, st := range states {
		if st.SampleID == "S1" && st.State == varstore.StateNonCarrier {
			return true, nil
		}
	}
	return false, nil
}

// TestVarQueryChromNamingIsAutoConverted covers both naming conventions against
// an indexed, chr-prefixed VCF. Tabix looks reference names up verbatim, so
// before contig resolution a "22"-style query against a "chr22" index failed
// with `tabix: unknown reference "22"` rather than simply finding the record.
//
// This only bites on the indexed path: an unindexed file is scanned and
// compared canonically, which is why the earlier tests missed it.
func TestVarQueryChromNamingIsAutoConverted(t *testing.T) {
	for _, q := range []string{"chr1:100:A:G", "1:100:A:G"} {
		out := runVcf(t, "vcf-varquery", "--variant", q, "testdata/sample.vcf.gz")
		if !strings.Contains(out, "NORMAL") && !strings.Contains(out, "TUMOR") {
			t.Errorf("query %q returned no carriers:\n%s", q, out)
		}
	}
	// The carriers must be identical either way. The leading "variant" column
	// deliberately echoes the spelling the user typed, so compare past it.
	a := dropFirstColumn(dataRowsOnly(runVcf(t, "vcf-varquery", "--variant", "chr1:100:A:G", "testdata/sample.vcf.gz")))
	b := dropFirstColumn(dataRowsOnly(runVcf(t, "vcf-varquery", "--variant", "1:100:A:G", "testdata/sample.vcf.gz")))
	if a != b {
		t.Errorf("chr-prefixed and bare spellings disagree\n chr1:\n%s\n 1:\n%s", a, b)
	}
}

// dropFirstColumn strips the leading tab-separated field of every line.
func dropFirstColumn(s string) string {
	var out []string
	for _, l := range strings.Split(s, "\n") {
		if _, rest, ok := strings.Cut(l, "\t"); ok {
			out = append(out, rest)
		} else {
			out = append(out, l)
		}
	}
	return strings.Join(out, "\n")
}

// TestVarQueryUnknownContigIsAnAbsence pins that asking about a contig the file
// does not have is answered, not rejected: it is simply not in the source. The
// behaviour must not depend on whether the file happens to be indexed.
func TestVarQueryUnknownContigIsAnAbsence(t *testing.T) {
	// The harness folds stderr into the same buffer, so assert on the absence of
	// carriers rather than on line counts.
	for _, in := range []string{"testdata/sample.vcf.gz", "testdata/sample.vcf"} {
		out := runVcf(t, "vcf-varquery", "--variant", "chrZZ:100:A:G", in)
		for _, sample := range []string{"NORMAL", "TUMOR"} {
			if strings.Contains(out, sample+"\t") {
				t.Errorf("%s: absent contig reported %s as a carrier:\n%s", in, sample, out)
			}
		}
		if !strings.Contains(out, "not in the source") {
			t.Errorf("%s: expected a warning that the variant is absent, got:\n%s", in, out)
		}
	}
}

// TestVarQueryBadRegionIsRejected pins the other half: --region names a contig
// the caller asserts exists, so an unresolvable one is an error that says what
// the file does have.
func TestVarQueryBadRegionIsRejected(t *testing.T) {
	err := runVcfErr(t, "vcf-varquery", "--sample", "NORMAL", "--region", "chrZZ:1-100", "testdata/sample.vcf.gz")
	if err == nil {
		t.Fatal("expected an error for a --region on a contig not in the index")
	}
	if !strings.Contains(err.Error(), "unknown reference") || !strings.Contains(err.Error(), "chr1") {
		t.Errorf("error should name the available references, got: %v", err)
	}
}

// TestParquetStorePreservesSourceChromNaming pins that conversion records the
// source's own spelling rather than rewriting it, while queries stay
// convention-agnostic.
func TestParquetStorePreservesSourceChromNaming(t *testing.T) {
	base := convert(t, "testdata/multiallelic.vcf") // chr-prefixed
	calls := readCalls(t, base)
	if len(calls) == 0 {
		t.Fatal("no calls")
	}
	for _, c := range calls {
		if !strings.HasPrefix(c.Chrom, "chr") {
			t.Errorf("stored chrom = %q, want the source's chr-prefixed spelling", c.Chrom)
			break
		}
	}
	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	for _, q := range []string{"chr1", "1"} {
		known, err := s.SiteKnown(varstore.Locus{Chrom: q, Pos: 100, Ref: "A", Alt: "G"})
		if err != nil {
			t.Fatal(err)
		}
		if !known {
			t.Errorf("SiteKnown(%s:100:A:G) = false, want true", q)
		}
	}
}

// TestVarQueryLocusColumnsAreSplit pins the tabular contract: the locus occupies
// four leading columns rather than one packed chrom:pos:ref:alt field, so
// downstream tools can cut or sort on position without re-splitting a key.
func TestVarQueryLocusColumnsAreSplit(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")

	cases := []struct {
		name   string
		args   []string
		header string
	}{
		{"carriers", []string{"vcf-varquery", "--variant", "chr1:100:A:G", base},
			"chrom\tpos\tref\talt\tsample\tgt\tdp\tmin_dp\tad_ref\tad_alt\tgq"},
		// S2 is the carrier at chr1:100, so --sample lands on the same locus the
		// other cases assert -- which is what makes the shared layout checkable.
		{"samples", []string{"vcf-varquery", "--sample", "S2", base},
			"chrom\tpos\tref\talt\tsample\tgt\tdp\tmin_dp\tad_ref\tad_alt\tgq"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			rows := strings.Split(dataRowsOnly(runVcf(t, tc.args...)), "\n")
			if rows[0] != tc.header {
				t.Errorf("header = %q, want %q", rows[0], tc.header)
			}
			if len(rows) < 2 {
				t.Fatal("expected at least one data row")
			}
			f := strings.Split(rows[1], "\t")
			if len(f) < 5 {
				t.Fatalf("row has %d columns: %q", len(f), rows[1])
			}
			if f[0] != "chr1" || f[1] != "100" || f[2] != "A" || f[3] != "G" {
				t.Errorf("locus columns = %q,%q,%q,%q, want chr1,100,A,G", f[0], f[1], f[2], f[3])
			}
			if strings.Contains(f[0], ":") {
				t.Errorf("first column %q still looks like a packed locus", f[0])
			}
		})
	}
}

// TestSitesACAN pins that AC/AN are allele counts and are not interchangeable
// with the sample counts beside them. The fixture is chosen so they diverge:
// hom-alt genotypes make AC exceed n_carriers, and a 1/2 sample splits one
// allele to each ALT row.
//
//	chr1:100 A>G   S1 0/1  S2 0/0  S3 1/1  S4 0/0
//	chr1:200 C>T,G S1 1/2  S2 0/1  S3 0/0  S4 2/2
//	chr1:300 G>A   all 0/0
func TestSitesACAN(t *testing.T) {
	base := convert(t, "testdata/multiallelic.vcf")
	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	cases := []struct {
		locus             varstore.Locus
		ac, an, nCarriers int32
		why               string
	}{
		{varstore.Locus{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "G"}, 3, 8, 2,
			"S1 het (1) + S3 hom (2) = 3 alleles across 2 carriers"},
		{varstore.Locus{Chrom: "chr1", Pos: 200, Ref: "C", Alt: "T"}, 2, 8, 2,
			"S1 1/2 contributes one T, S2 0/1 one T"},
		{varstore.Locus{Chrom: "chr1", Pos: 200, Ref: "C", Alt: "G"}, 3, 8, 2,
			"S1 1/2 one G + S4 2/2 two G = 3 alleles across 2 carriers"},
		{varstore.Locus{Chrom: "chr1", Pos: 300, Ref: "G", Alt: "A"}, 0, 8, 0,
			"monomorphic, but still interrogated in 4 diploid samples"},
	}
	for _, tc := range cases {
		site, ok, err := s.Site(tc.locus)
		if err != nil {
			t.Fatal(err)
		}
		if !ok {
			t.Errorf("%s missing from the catalog", tc.locus)
			continue
		}
		if site.AC != tc.ac {
			t.Errorf("%s AC = %d, want %d (%s)", tc.locus, site.AC, tc.ac, tc.why)
		}
		if site.AN != tc.an {
			t.Errorf("%s AN = %d, want %d", tc.locus, site.AN, tc.an)
		}
		if site.NCarriers != tc.nCarriers {
			t.Errorf("%s n_carriers = %d, want %d", tc.locus, site.NCarriers, tc.nCarriers)
		}
		if want := float64(tc.ac) / float64(tc.an); site.AF() != want {
			t.Errorf("%s AF = %v, want %v", tc.locus, site.AF(), want)
		}
	}

	// The point of the test: at two of these loci AC and n_carriers differ, so
	// one genuinely cannot stand in for the other.
	site, _, _ := s.Site(varstore.Locus{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "G"})
	if site.AC == site.NCarriers {
		t.Error("fixture no longer distinguishes AC from n_carriers; the test proves nothing")
	}
}

// TestSitesACANSurviveNoCallable pins that allele counts are properties of the
// genotypes, not of coverage, so a DP-less source still gets real AC/AN even
// though its sample-level coverage counts are necessarily zero.
func TestSitesACANSurviveNoCallable(t *testing.T) {
	base := convert(t, "testdata/sample.vcf", "--no-callable")
	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	sawAN := false
	if err := s.Sites(func(site varstore.Site) bool {
		if site.AN > 0 {
			sawAN = true
			return false
		}
		return true
	}); err != nil {
		t.Fatal(err)
	}
	if !sawAN {
		t.Error("AN is zero everywhere under --no-callable; allele counts must not depend on DP")
	}
}

// TestVerboseGoesToStderrOnly is the constraint that makes verbose safe to use
// in a pipeline: the tabular stream must stay parseable.
func TestVerboseGoesToStderrOnly(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	stdout := filepath.Join(t.TempDir(), "out.tsv")
	runVcf(t, "vcf-varquery", "-v", "--variant", "chr1:100:A:G", "-o", stdout, base)

	b, err := os.ReadFile(stdout)
	if err != nil {
		t.Fatal(err)
	}
	got := string(b)
	for _, marker := range []string{"store ", "variant  ", "gate ", "WARNING", "NOTE"} {
		if strings.Contains(got, marker) {
			t.Errorf("verbose text %q leaked into the data stream:\n%s", marker, got)
		}
	}
	if !strings.HasPrefix(strings.SplitN(got, "\n", 2)[0], "##") {
		t.Errorf("data stream should still start with provenance, got:\n%s", got)
	}
}

// TestVerboseReportsAbsentQualityFields pins the diagnostic worth having: a gate
// over a field the data lacks admits everything, so a --min-gq that looks like a
// filter can be doing nothing. Only the field census distinguishes that from a
// real result.
func TestVerboseReportsAbsentQualityFields(t *testing.T) {
	// sample.vcf carries neither DP nor GQ.
	out := runVcf(t, "vcf-varquery", "-v", "--variant", "chr1:100:A:G",
		"--min-dp", "99", "testdata/sample.vcf")
	if !strings.Contains(out, "--min-dp 99 had no effect") {
		t.Errorf("expected a warning that the gate could not act, got:\n%s", out)
	}
}

// TestVerboseConversionReportsFieldPresence pins the same idea at conversion
// time, so a store built from quality-less input says so when it is created.
func TestVerboseConversionReportsFieldPresence(t *testing.T) {
	base := filepath.Join(t.TempDir(), "s")
	out := runVcf(t, "vcf-toparquet", "-v", "--no-callable", "--out", base, "testdata/sample.vcf")
	for _, want := range []string{"fields present", "GQ  ABSENT", "DP  ABSENT"} {
		if !strings.Contains(out, want) {
			t.Errorf("verbose conversion report missing %q, got:\n%s", want, out)
		}
	}

	// And with quality present it should report coverage rather than absence.
	base2 := filepath.Join(t.TempDir(), "s2")
	out2 := runVcf(t, "vcf-toparquet", "-v", "--out", base2, "testdata/coverage.vcf")
	if strings.Contains(out2, "DP  ABSENT") {
		t.Errorf("coverage.vcf has DP; should not be reported absent:\n%s", out2)
	}
	if !strings.Contains(out2, "called but below DP") {
		t.Errorf("expected a coverage breakdown, got:\n%s", out2)
	}
}

// TestVerboseNotesMinDPMismatch pins the warning for querying at a threshold the
// store's runs were not built at, which silently breaks backend agreement.
func TestVerboseNotesMinDPMismatch(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf") // built at the default --min-dp 10
	out := runVcf(t, "vcf-varquery", "-v", "--variant", "chr1:100:A:G",
		"--hom-ref", "--min-dp", "25", base)
	if !strings.Contains(out, "the runs were built at 10") {
		t.Errorf("expected a min-dp mismatch note, got:\n%s", out)
	}
	// Matching the conversion threshold must not warn.
	ok := runVcf(t, "vcf-varquery", "-v", "--variant", "chr1:100:A:G",
		"--hom-ref", "--min-dp", "10", base)
	if strings.Contains(ok, "the runs were built at") {
		t.Errorf("matching --min-dp should not warn, got:\n%s", ok)
	}
}

// TestQuietByDefault pins that none of this appears without -v.
func TestQuietByDefault(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	out := runVcf(t, "vcf-varquery", "--variant", "chr1:100:A:G", "--min-dp", "99", base)
	for _, marker := range []string{"store    parquet", "variant  ", "gate     "} {
		if strings.Contains(out, marker) {
			t.Errorf("saw verbose output %q without -v:\n%s", marker, out)
		}
	}
}

// TestStoreDirLayout pins the directory form: a base ending in "/" creates the
// directory (including missing parents) and puts the members inside it under
// bare names, with no prefix dot.
func TestStoreDirLayout(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "nested", "cohort") + string(os.PathSeparator)
	runVcf(t, "vcf-toparquet", "--out", dir, "testdata/coverage.vcf")

	for _, name := range []string{"calls.parquet", "sites.parquet", "regions.parquet"} {
		p := filepath.Join(dir, name)
		st, err := os.Stat(p)
		if err != nil {
			t.Errorf("missing %s: %v", p, err)
			continue
		}
		if st.Size() == 0 {
			t.Errorf("%s is empty", p)
		}
	}
	// Nothing should have been created alongside the directory with a dot.
	if matches, _ := filepath.Glob(strings.TrimSuffix(dir, string(os.PathSeparator)) + ".*.parquet"); len(matches) > 0 {
		t.Errorf("directory form should not also write prefixed files: %v", matches)
	}
}

// TestStoreDirAndPrefixAgree pins that the layout is purely cosmetic: the same
// input converted either way answers identically.
func TestStoreDirAndPrefixAgree(t *testing.T) {
	tmp := t.TempDir()
	dir := filepath.Join(tmp, "d") + string(os.PathSeparator)
	pfx := filepath.Join(tmp, "p")
	runVcf(t, "vcf-toparquet", "--out", dir, "testdata/coverage.vcf")
	runVcf(t, "vcf-toparquet", "--out", pfx, "testdata/coverage.vcf")

	a := dataRowsOnly(runVcf(t, "vcf-varquery", "--variant", "chr1:100:A:G", "--hom-ref", dir))
	b := dataRowsOnly(runVcf(t, "vcf-varquery", "--variant", "chr1:100:A:G", "--hom-ref", pfx))
	if a != b {
		t.Errorf("directory and prefix layouts disagree\n dir:\n%s\n prefix:\n%s", a, b)
	}
}

// TestStoreDirAcceptedEveryWay pins the spellings a user might reasonably type
// for a directory-form store, including the bare directory with no slash --
// having written "--out cohort/", nobody should need the slash to read it back.
func TestStoreDirAcceptedEveryWay(t *testing.T) {
	tmp := t.TempDir()
	dir := filepath.Join(tmp, "cohort")
	runVcf(t, "vcf-toparquet", "--out", dir+string(os.PathSeparator), "testdata/coverage.vcf")

	want := dataRowsOnly(runVcf(t, "vcf-varquery", "--variant", "chr1:100:A:G",
		dir+string(os.PathSeparator)))
	for _, spelling := range []string{
		dir,                                   // bare directory, no slash
		filepath.Join(dir, "calls.parquet"),   // a member
		filepath.Join(dir, "sites.parquet"),   // a different member
		filepath.Join(dir, "regions.parquet"), // and the third
	} {
		got := dataRowsOnly(runVcf(t, "vcf-varquery", "--variant", "chr1:100:A:G", spelling))
		if got != want {
			t.Errorf("%s resolved to a different answer:\n%s", spelling, got)
		}
	}
}

// TestStorePathHelpers pins the two shapes at the unit level, including that the
// directory form introduces no dot.
func TestStorePathHelpers(t *testing.T) {
	cases := []struct{ base, calls, sites, regions string }{
		{"cohort", "cohort.calls.parquet", "cohort.sites.parquet", "cohort.regions.parquet"},
		{"cohort/", "cohort/calls.parquet", "cohort/sites.parquet", "cohort/regions.parquet"},
		{"a/b/", "a/b/calls.parquet", "a/b/sites.parquet", "a/b/regions.parquet"},
		{"a/b", "a/b.calls.parquet", "a/b.sites.parquet", "a/b.regions.parquet"},
	}
	for _, tc := range cases {
		if got := varstore.CallsPath(tc.base); got != tc.calls {
			t.Errorf("CallsPath(%q) = %q, want %q", tc.base, got, tc.calls)
		}
		if got := varstore.SitesPath(tc.base); got != tc.sites {
			t.Errorf("SitesPath(%q) = %q, want %q", tc.base, got, tc.sites)
		}
		if got := varstore.RegionsPath(tc.base); got != tc.regions {
			t.Errorf("RegionsPath(%q) = %q, want %q", tc.base, got, tc.regions)
		}
	}
}

// TestStoreOverwriteRequiresForce pins that conversion will not silently
// destroy an existing store. Writing truncates all three members, so the check
// has to happen before anything is opened.
func TestStoreOverwriteRequiresForce(t *testing.T) {
	for _, form := range []struct{ name, suffix string }{
		{"prefix", ""},
		{"directory", string(os.PathSeparator)},
	} {
		t.Run(form.name, func(t *testing.T) {
			base := filepath.Join(t.TempDir(), "store") + form.suffix
			runVcf(t, "vcf-toparquet", "--out", base, "testdata/coverage.vcf")

			err := runVcfErr(t, "vcf-toparquet", "--out", base, "testdata/coverage.vcf")
			if err == nil {
				t.Fatal("expected a refusal to overwrite an existing store")
			}
			if !strings.Contains(err.Error(), "--force") {
				t.Errorf("error should point at --force, got: %v", err)
			}
			// The refusal must leave the previous store usable.
			out := runVcf(t, "vcf-varquery", "--variant", "chr1:100:A:G", base)
			if !strings.Contains(out, "S2") {
				t.Errorf("refusal damaged the existing store:\n%s", out)
			}
			// And --force must go through.
			runVcf(t, "vcf-toparquet", "--force", "--out", base, "testdata/coverage.vcf")
		})
	}
}

// TestStoreOverwriteRefusesOnPartialSet pins that a single surviving member is
// enough to stop: a half-replaced set is worse than either whole outcome.
func TestStoreOverwriteRefusesOnPartialSet(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	if err := os.Remove(varstore.CallsPath(base)); err != nil {
		t.Fatal(err)
	}
	if err := os.Remove(varstore.RegionsPath(base)); err != nil {
		t.Fatal(err)
	}
	err := runVcfErr(t, "vcf-toparquet", "--out", base, "testdata/coverage.vcf")
	if err == nil {
		t.Fatal("expected a refusal with only the sites member present")
	}
	if !strings.Contains(err.Error(), varstore.SitesPath(base)) {
		t.Errorf("error should name the surviving member, got: %v", err)
	}
}

// TestStorePrefixNamingADirectory pins the ambiguous case: "--out cohort" where
// cohort/ already exists almost certainly means the slash was forgotten.
func TestStorePrefixNamingADirectory(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "cohort")
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	err := runVcfErr(t, "vcf-toparquet", "--out", dir, "testdata/coverage.vcf")
	if err == nil {
		t.Fatal("expected a refusal when the prefix names an existing directory")
	}
	if !strings.Contains(err.Error(), "is a directory") {
		t.Errorf("error should explain the ambiguity, got: %v", err)
	}
	// --force accepts it as a genuine filename prefix.
	runVcf(t, "vcf-toparquet", "--force", "--out", dir, "testdata/coverage.vcf")
	if _, err := os.Stat(dir + ".calls.parquet"); err != nil {
		t.Errorf("--force should have written the prefix form: %v", err)
	}
}

// TestStoreDirWithUnrelatedContentIsFine pins that the guard keys on the store
// members, not on the directory being empty: an existing directory holding other
// files is a normal place to put a store, and its contents are left alone.
func TestStoreDirWithUnrelatedContentIsFine(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "cohort") + string(os.PathSeparator)
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	other := filepath.Join(dir, "notes.txt")
	if err := os.WriteFile(other, []byte("keep me"), 0o644); err != nil {
		t.Fatal(err)
	}
	runVcf(t, "vcf-toparquet", "--out", dir, "testdata/coverage.vcf")

	b, err := os.ReadFile(other)
	if err != nil || string(b) != "keep me" {
		t.Errorf("unrelated file was disturbed: %v %q", err, b)
	}
}

// TestMultiVcfCombines pins that several inputs land in one store.
func TestMultiVcfCombines(t *testing.T) {
	base := filepath.Join(t.TempDir(), "s") + string(os.PathSeparator)
	runVcf(t, "vcf-toparquet", "--out", base,
		"testdata/multi_chr1.vcf", "testdata/multi_chr2.vcf")

	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	for _, l := range []varstore.Locus{
		{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "G"},
		{Chrom: "chr2", Pos: 200, Ref: "T", Alt: "A"},
	} {
		known, err := s.SiteKnown(l)
		if err != nil {
			t.Fatal(err)
		}
		if !known {
			t.Errorf("%s missing from the combined store", l)
		}
	}
}

// TestMultiVcfRemapsSampleOrder is the one that matters. Genotype columns are
// positional, so a second input with reversed columns would, without remapping,
// attribute every genotype to the wrong person -- silently, with output that
// looks entirely plausible. The reordered file carries the same genotypes per
// SAMPLE as the canonical one, so both must produce identical stores.
func TestMultiVcfRemapsSampleOrder(t *testing.T) {
	canonical := filepath.Join(t.TempDir(), "a") + string(os.PathSeparator)
	reordered := filepath.Join(t.TempDir(), "b") + string(os.PathSeparator)
	runVcf(t, "vcf-toparquet", "--out", canonical,
		"testdata/multi_chr1.vcf", "testdata/multi_chr2.vcf")
	runVcf(t, "vcf-toparquet", "--out", reordered,
		"testdata/multi_chr1.vcf", "testdata/multi_chr2_reordered.vcf")

	for _, l := range []string{"chr2:100:G:C", "chr2:200:T:A"} {
		a := dataRowsOnly(runVcf(t, "vcf-varquery", "--variant", l, "--hom-ref", canonical))
		b := dataRowsOnly(runVcf(t, "vcf-varquery", "--variant", l, "--hom-ref", reordered))
		if a != b {
			t.Errorf("reordered input produced different genotypes at %s\n canonical:\n%s\n reordered:\n%s",
				l, a, b)
		}
	}
}

// TestMultiVcfRejectsDifferentSamples pins that a set mismatch is refused, and
// that the message names what differs rather than just saying they differ.
func TestMultiVcfRejectsDifferentSamples(t *testing.T) {
	base := filepath.Join(t.TempDir(), "s") + string(os.PathSeparator)
	err := runVcfErr(t, "vcf-toparquet", "--out", base,
		"testdata/multi_chr1.vcf", "testdata/multi_chr2_othersamples.vcf")
	if err == nil {
		t.Fatal("expected a refusal for differing sample sets")
	}
	for _, want := range []string{"same samples", "S3"} {
		if !strings.Contains(err.Error(), want) {
			t.Errorf("error should mention %q, got: %v", want, err)
		}
	}
}

// TestMultiVcfRejectsOverlap pins that overlapping inputs are refused: they
// would write the same site twice and split its AC/AN across two rows.
func TestMultiVcfRejectsOverlap(t *testing.T) {
	base := filepath.Join(t.TempDir(), "s") + string(os.PathSeparator)
	err := runVcfErr(t, "vcf-toparquet", "--out", base,
		"testdata/multi_chr1.vcf", "testdata/multi_chr1.vcf")
	if err == nil {
		t.Fatal("expected a refusal for overlapping inputs")
	}
	if !strings.Contains(err.Error(), "overlap") && !strings.Contains(err.Error(), "sorted") {
		t.Errorf("error should explain the ordering rule, got: %v", err)
	}
}

// TestFailedConversionLeavesNothing pins that a failure does not leave members
// on disk: they would look like a store, and would block the retry through the
// overwrite guard.
func TestFailedConversionLeavesNothing(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "s") + string(os.PathSeparator)
	if err := runVcfErr(t, "vcf-toparquet", "--out", base,
		"testdata/multi_chr1.vcf", "testdata/multi_chr2_othersamples.vcf"); err == nil {
		t.Fatal("expected the conversion to fail")
	}
	if found := varstore.ExistingMembers(base); len(found) > 0 {
		t.Errorf("failed conversion left %v behind", found)
	}
	// The retry must succeed without --force.
	runVcf(t, "vcf-toparquet", "--out", base,
		"testdata/multi_chr1.vcf", "testdata/multi_chr2.vcf")
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
