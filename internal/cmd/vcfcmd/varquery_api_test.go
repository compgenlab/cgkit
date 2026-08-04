package vcfcmd

import (
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/varstore"
)

// Tests for the varstore.Calls(Query) surface, exercised through both backends.
//
// The point of the single-method reshape is that site selection and sample
// selection are independent axes, so these cover the combinations the four old
// methods could not express -- several loci at once, several samples at once, and
// both axes together. They live here rather than in varstore because the fixtures
// and the VCF-to-store converter do.

// bothBackends runs a query against a VCF and against a store converted from it,
// asserts they agree, and returns the rows as identity strings.
//
// Quality columns are deliberately excluded from the comparison: a
// store-reconstructed reference call has no DP/AD/GQ where a VCF does. What must
// match is which rows appear, and in what order.
func bothBackends(t *testing.T, vcfPath string, q varstore.Query) []string {
	t.Helper()
	base := convert(t, vcfPath)

	v, err := varstore.OpenVcf(vcfPath)
	if err != nil {
		t.Fatalf("OpenVcf: %v", err)
	}
	defer v.Close()
	p, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatalf("OpenParquet: %v", err)
	}
	defer p.Close()

	keys := func(s varstore.Store, which string) []string {
		calls, err := varstore.CollectCalls(s, q)
		if err != nil {
			t.Fatalf("%s Calls: %v", which, err)
		}
		out := make([]string, 0, len(calls))
		for _, c := range calls {
			out = append(out, c.SampleID+" "+c.Locus().String()+" "+c.GT)
		}
		return out
	}
	a, b := keys(v, "vcf"), keys(p, "parquet")
	if strings.Join(a, "|") != strings.Join(b, "|") {
		t.Errorf("backends disagree\n vcf:     %v\n parquet: %v", a, b)
	}
	return b
}

// TestQuerySeveralLociAtOnce is the case the old API could only serve by looping,
// which is what made a panel query cost one row-group decode per variant.
func TestQuerySeveralLociAtOnce(t *testing.T) {
	got := bothBackends(t, "testdata/coverage.vcf", varstore.Query{Loci: []varstore.Locus{
		{Chrom: "chr1", Pos: 500, Ref: "A", Alt: "T"},
		{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "G"},
	}})
	// Store order, not the order the loci were named.
	want := "S2 chr1:100:A:G 0/1|S1 chr1:500:A:T 0/1"
	if strings.Join(got, "|") != want {
		t.Errorf("got %v\nwant %s (store order, not query order)", got, want)
	}
}

// TestQueryBothAxes pins the combination that used to be a hard error at the CLI
// and inexpressible in the library: these samples, at these sites.
func TestQueryBothAxes(t *testing.T) {
	locus := []varstore.Locus{{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "G"}}

	got := bothBackends(t, "testdata/coverage.vcf", varstore.Query{
		Loci:    locus,
		Samples: []string{"S1"},
	})
	if len(got) != 0 {
		t.Errorf("S1 does not carry chr1:100 A>G, so an ALT-only query is empty; got %v", got)
	}

	// With reference calls the same query answers "what is S1's genotype there".
	got = bothBackends(t, "testdata/coverage.vcf", varstore.Query{
		Loci:       locus,
		Samples:    []string{"S1"},
		Gate:       varstore.Gate{MinDP: 10},
		IncludeRef: true,
	})
	if len(got) != 1 || !strings.HasPrefix(got[0], "S1 chr1:100:A:G 0/0") {
		t.Errorf("want S1 reported 0/0 at chr1:100, got %v", got)
	}
}

// TestQueryEmptySelectorMeansEverything pins that an empty axis is "no
// restriction" rather than "nothing", so both axes read the same way.
func TestQueryEmptySelectorMeansEverything(t *testing.T) {
	all := bothBackends(t, "testdata/coverage.vcf", varstore.Query{})
	if len(all) != 2 {
		t.Errorf("the zero Query should return every ALT call (2 in this fixture), got %d: %v",
			len(all), all)
	}
	n := 0
	for _, s := range []string{"S1", "S2"} {
		n += len(bothBackends(t, "testdata/coverage.vcf", varstore.Query{Samples: []string{s}}))
	}
	if n != len(all) {
		t.Errorf("naming each sample gave %d rows, the unrestricted query %d", n, len(all))
	}
}

// TestQuerySiblingAltDisqualifiesReference is the regression test for a bug the
// reshape introduced and the existing suite caught.
//
// At chr1:200 C>T,G sample S2 is 0/1 -- it carries T. Asking only about the G
// locus, the evidence that S2 is not reference lives at the SIBLING locus, so a
// site-level exclusion test loses it and fabricates a 0/0. The exclusion has to be
// keyed on the source record, which is what plan.wantsRecord is for.
func TestQuerySiblingAltDisqualifiesReference(t *testing.T) {
	got := bothBackends(t, "testdata/multiallelic.vcf", varstore.Query{
		Loci:       []varstore.Locus{{Chrom: "chr1", Pos: 200, Ref: "C", Alt: "G"}},
		IncludeRef: true,
	})
	refs := 0
	for _, row := range got {
		if !strings.HasSuffix(row, varstore.HomRefGT) {
			continue
		}
		refs++
		if strings.HasPrefix(row, "S2 ") {
			t.Errorf("S2 carries T at that record, so it is not reference at the G locus: %q\nrows: %v",
				row, got)
		}
	}
	if refs != 1 {
		t.Errorf("only S3 is all-reference at chr1:200; got %d reference rows: %v", refs, got)
	}
}

// TestQueryIncludeRefRefusesAtSetup pins that an unclassifiable store fails when
// Calls is called, before any row is yielded, so a caller cannot mistake a
// silently ALT-only stream for a complete answer.
func TestQueryIncludeRefRefusesAtSetup(t *testing.T) {
	// --no-callable rather than a deleted regions file: since the manifest
	// records what each member held, removing one is corruption and is caught
	// at open. This is about the store that legitimately tracked no coverage.
	base := filepath.Join(t.TempDir(), "store")
	runVcf(t, "vcf-toparquet", "--no-callable", "--out", base, "testdata/coverage.vcf")
	p, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer p.Close()

	if _, err := p.Calls(varstore.Query{IncludeRef: true}); err == nil {
		t.Error("expected ErrNotClassifiable from Calls itself, not a silent ALT-only stream")
	}
	if _, err := p.Calls(varstore.Query{}); err != nil {
		t.Errorf("an ALT-only query needs no regions file: %v", err)
	}
}

// TestQueryStopsEarly pins that abandoning the iterator stops the scan, which is
// what the streaming form exists for.
func TestQueryStopsEarly(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")
	p, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer p.Close()

	seq, err := p.Calls(varstore.Query{})
	if err != nil {
		t.Fatal(err)
	}
	n := 0
	for _, err := range seq {
		if err != nil {
			t.Fatal(err)
		}
		n++
		break
	}
	if n != 1 {
		t.Errorf("breaking out of the iterator should stop after one row, got %d", n)
	}
}

// TestQueryRowsAreSampleOrderedWithinALocus pins the ordering contract, which the
// first streaming implementation broke.
//
// Rows arrive ordered by (chrom, pos, ref, alt, sample). Emitting all the ALT
// calls for a locus and then all its reference calls satisfies the first four but
// not the fifth: at contigs.vcf chr10:200 S1 is 0/0 and S2 is 0/1, so the ALT
// block put S2 ahead of S1.
func TestQueryRowsAreSampleOrderedWithinALocus(t *testing.T) {
	got := bothBackends(t, "testdata/contigs.vcf", varstore.Query{
		Loci:       []varstore.Locus{{Chrom: "chr10", Pos: 200, Ref: "T", Alt: "C"}},
		Gate:       varstore.Gate{MinDP: 10},
		IncludeRef: true,
	})
	if len(got) != 2 {
		t.Fatalf("want a row per sample, got %v", got)
	}
	if !strings.HasPrefix(got[0], "S1 ") || !strings.HasPrefix(got[1], "S2 ") {
		t.Errorf("rows not in sample order: %v\nS1 is the reference call and S2 the ALT "+
			"call, so a type-major order would put S2 first", got)
	}
}
