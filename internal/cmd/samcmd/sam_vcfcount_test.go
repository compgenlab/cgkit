package samcmd

import (
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio"
	"github.com/compgenlab/cghts/vcf"
	"github.com/spf13/pflag"
)

// Reference used by the synthetic reads below, 1-based positions 1..20:
//
//	pos  1234567890123456789 0
//	ref  ACGTACGTACGTACGTACGT
const vcRef = "ACGTACGTACGTACGTACGT"

func vcVariantFor(t *testing.T, pos int, ref string, alts ...string) *vcVariant {
	t.Helper()
	line := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t" + strconv.Itoa(pos) + "\t.\t" + ref + "\t" + strings.Join(alts, ",") + "\t.\tPASS\t.\n"
	r, err := vcf.NewVcfReader(strings.NewReader("##fileformat=VCFv4.2\n" + line))
	if err != nil {
		t.Fatal(err)
	}
	rec, err := r.NextRecord()
	if err != nil {
		t.Fatal(err)
	}
	v := newVcVariant(rec)
	if v == nil {
		t.Fatalf("variant %d %s>%v not countable", pos, ref, alts)
	}
	v.bamRef = "chr1"
	return v
}

// read builds an aligned read at 1-based pos with the given cigar and seq.
// qual may be "" for no qualities.
func read(name string, flag, pos int, cigar, seq, qual string) *alignedRead {
	if qual == "" {
		qual = "*"
	}
	return newAlignedRead(&htsio.SamRecord{ReadName: name, Flag: flag, RefName: "chr1", Pos: pos, Cigar: cigar, Seq: seq, Qual: qual})
}

func count(t *testing.T, v *vcVariant, reads []*alignedRead, minBQ int, once bool) (ref int, alts []int, total int) {
	t.Helper()
	countVariants([]*vcVariant{v}, reads, minBQ, once)
	return v.refCount, v.altCounts, v.total
}

func TestVcfCountSNV(t *testing.T) {
	v := vcVariantFor(t, 6, "C", "T")
	reads := []*alignedRead{
		read("r1", 0, 1, "10M", "ACGTACGTAC", ""), // C at 6: ref
		read("r2", 0, 2, "10M", "CGTATGTACG", ""), // T at 6: alt
		read("r3", 0, 3, "10M", "GTAGGTACGT", ""), // G at 6: other
		read("r5", 0, 1, "5M", "ACGTA", ""),       // ends before
		read("r6", 0, 1, "5M2D5M", "ACGTAGTACG", ""),
		read("r4", 0, 7, "10M", "GTACGTACGT", ""), // starts after: not spanning
	}
	ref, alts, total := count(t, v, reads, 13, true)
	if ref != 1 || alts[0] != 1 || total != 4 {
		t.Errorf("ref=%d alt=%d total=%d, want 1 1 4 (r6 deletes the site: other)", ref, alts[0], total)
	}
}

func TestVcfCountBaseQuality(t *testing.T) {
	v := vcVariantFor(t, 6, "C", "T")
	reads := []*alignedRead{
		read("hi", 0, 1, "10M", "ACGTATGTAC", "IIIIIIIIII"),
		read("lo", 0, 1, "10M", "ACGTATGTAC", "IIIII#IIII"), // Q2 at the site
	}
	_, alts, total := count(t, v, reads, 13, true)
	if alts[0] != 1 || total != 1 {
		t.Errorf("alt=%d total=%d, want 1 1: a low-quality base is not counted", alts[0], total)
	}
}

func TestVcfCountDeletion(t *testing.T) {
	// ref positions 4..6 = T A C; deletion of AC after anchor T.
	v := vcVariantFor(t, 4, "TAC", "T")
	reads := []*alignedRead{
		read("ref", 0, 1, "10M", "ACGTACGTAC", ""),
		read("alt", 0, 1, "4M2D6M", "ACGTGTACGT", ""),
		read("longer", 0, 1, "4M3D6M", "ACGTTACGTA", ""),  // a different deletion: other
		read("endsOn", 0, 1, "6M", "ACGTAC", ""),          // no base after REF: not counted
		read("anchorDel", 0, 1, "3M3D5M", "ACGCGTAC", ""), // anchor deleted: other
	}
	ref, alts, total := count(t, v, reads, 13, true)
	if ref != 1 || alts[0] != 1 || total != 4 {
		t.Errorf("ref=%d alt=%d total=%d, want 1 1 4", ref, alts[0], total)
	}
}

func TestVcfCountInsertion(t *testing.T) {
	v := vcVariantFor(t, 4, "T", "TGG")
	reads := []*alignedRead{
		read("ref", 0, 1, "8M", "ACGTACGT", ""),
		read("alt", 0, 1, "4M2I4M", "ACGTGGACGT", ""),
		read("otherIns", 0, 1, "4M1I4M", "ACGTGACGT", ""),
		read("endsOnAnchor", 0, 1, "4M", "ACGT", ""), // cannot tell
	}
	ref, alts, total := count(t, v, reads, 13, true)
	if ref != 1 || alts[0] != 1 || total != 3 {
		t.Errorf("ref=%d alt=%d total=%d, want 1 1 3", ref, alts[0], total)
	}
}

func TestVcfCountMultiallelic(t *testing.T) {
	v := vcVariantFor(t, 6, "C", "T", "G", "<NON_REF>")
	reads := []*alignedRead{
		read("a", 0, 1, "10M", "ACGTACGTAC", ""),
		read("b", 0, 1, "10M", "ACGTATGTAC", ""),
		read("c", 0, 1, "10M", "ACGTAGGTAC", ""),
		read("d", 0, 1, "10M", "ACGTAGGTAC", ""),
	}
	ref, alts, total := count(t, v, reads, 13, true)
	if ref != 1 || alts[0] != 1 || alts[1] != 2 || alts[2] != 0 || total != 4 {
		t.Errorf("ref=%d alts=%v total=%d, want 1 [1 2 0] 4", ref, alts, total)
	}
	if v.alts[2] != "" {
		t.Errorf("<NON_REF> should be uncountable")
	}
}

func TestVcfCountMates(t *testing.T) {
	v := vcVariantFor(t, 6, "C", "T")
	reads := []*alignedRead{
		read("agree", 0x1, 1, "10M", "ACGTATGTAC", ""),
		read("agree", 0x1, 2, "10M", "CGTATGTACG", ""),
		read("conflict", 0x1, 1, "10M", "ACGTACGTAC", ""),
		read("conflict", 0x1, 3, "10M", "GTATGTACGT", ""),
		read("single", 0x1, 1, "10M", "ACGTACGTAC", ""),
	}
	ref, alts, total := count(t, v, reads, 13, true)
	if ref != 1 || alts[0] != 1 || total != 2 {
		t.Errorf("once per fragment: ref=%d alt=%d total=%d, want 1 1 2", ref, alts[0], total)
	}
	v = vcVariantFor(t, 6, "C", "T")
	ref, alts, total = count(t, v, reads, 13, false)
	if ref != 2 || alts[0] != 3 || total != 5 {
		t.Errorf("both mates: ref=%d alt=%d total=%d, want 2 3 5", ref, alts[0], total)
	}
}

func TestVcfCountSplicedAcross(t *testing.T) {
	v := vcVariantFor(t, 6, "C", "T")
	reads := []*alignedRead{read("s", 0, 1, "4M5N4M", "ACGTCGTA", "")}
	_, _, total := count(t, v, reads, 13, true)
	if total != 0 {
		t.Errorf("total=%d, want 0: a read spliced across the site does not cover it", total)
	}
}

func TestIsHet(t *testing.T) {
	for gt, want := range map[string]bool{"0/1": true, "1|0": true, "1/2": true, "0/0": false, "1/1": false, "./1": false, "0": false} {
		if got := isHet(gt); got != want {
			t.Errorf("isHet(%q)=%v, want %v", gt, got, want)
		}
	}
}

// End to end against testdata/test.bam: read1 is chr1:100-149
// ACGTACGTAC repeated, so chr1:101 is C and chr1:103 is T; read5 is chr2:200-249.
func runVcfCount(t *testing.T, vcfBody string, args ...string) (string, error) {
	t.Helper()
	dir := t.TempDir()
	vcfPath := filepath.Join(dir, "in.vcf")
	hdr := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
	if err := os.WriteFile(vcfPath, []byte(hdr+vcfBody), 0o644); err != nil {
		t.Fatal(err)
	}
	samVcfCountCmd.Flags().VisitAll(func(f *pflag.Flag) {
		if f.Changed {
			f.Value.Set(f.DefValue)
			f.Changed = false
		}
	})
	var out, errb bytes.Buffer
	samVcfCountCmd.SetOut(&out)
	samVcfCountCmd.SetErr(&errb)
	if err := samVcfCountCmd.ParseFlags(args); err != nil {
		t.Fatal(err)
	}
	err := samVcfCountCmd.RunE(samVcfCountCmd, []string{"testdata/test.bam", vcfPath})
	return dataRows(out.String()), err
}

func dataRows(s string) string {
	var rows []string
	for _, l := range strings.Split(strings.TrimSpace(s), "\n") {
		if !strings.HasPrefix(l, "##") {
			rows = append(rows, l)
		}
	}
	return strings.Join(rows, "\n")
}

func TestSamVcfCountEndToEnd(t *testing.T) {
	body := "" +
		"chr1\t101\trs1\tC\tT\t.\tPASS\t.\tGT\t0/1\n" +
		"1\t103\t.\tG\tT\t.\tPASS\t.\tGT\t1/1\n" + // Ensembl name; the read carries T = alt
		"chr1\t104\t.\tT\t<DEL>\t.\tPASS\t.\tGT\t0/1\n" + // symbolic: skipped
		"chr1\t300\t.\tA\tG\t.\tlowq\t.\tGT\t0/1\n" + // no reads
		"chr2\t210\t.\tA\tC,G\t.\tPASS\t.\tGT\t1/2\n"
	got, err := runVcfCount(t, body, "--id", "--total", "--gt")
	if err != nil {
		t.Fatal(err)
	}
	want := "chrom\tpos\tid\tref\talt\tref_count\talt_count\ttotal_count\tgt\n" +
		"chr1\t101\trs1\tC\tT\t1\t0\t1\t0/1\n" +
		"1\t103\t.\tG\tT\t0\t1\t1\t1/1\n" +
		"chr1\t300\t.\tA\tG\t0\t0\t0\t0/1\n" +
		"chr2\t210\t.\tA\tC\t1\t0\t1\t1/2\n" +
		"chr2\t210\t.\tA\tG\t1\t0\t1\t1/2"
	if got != want {
		t.Errorf("got\n%s\nwant\n%s", got, want)
	}

	got, err = runVcfCount(t, body, "--passing", "--het")
	if err != nil {
		t.Fatal(err)
	}
	if n := strings.Count(got, "\n"); n != 3 { // header + 101 + two chr2 rows
		t.Errorf("--passing --het kept %d rows, want 3:\n%s", n, got)
	}
}

func TestSamVcfCountMissingContig(t *testing.T) {
	body := "chr3\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\nchr1\t101\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n"
	if _, err := runVcfCount(t, body); err == nil || !strings.Contains(err.Error(), "--skip-missing") {
		t.Errorf("err = %v, want a missing-contig error naming --skip-missing", err)
	}
	got, err := runVcfCount(t, body, "--skip-missing")
	if err != nil {
		t.Fatal(err)
	}
	if strings.Contains(got, "chr3") || !strings.Contains(got, "chr1\t101") {
		t.Errorf("--skip-missing output:\n%s", got)
	}
}

func TestSamVcfCountMapQFilter(t *testing.T) {
	body := "chr1\t101\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n"
	got, err := runVcfCount(t, body, "--min-mapq", "61")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.HasSuffix(got, "chr1\t101\tC\tT\t0\t0") {
		t.Errorf("--min-mapq 61 should drop the MAPQ-60 read:\n%s", got)
	}
}

func TestSamVcfCountMinDepth(t *testing.T) {
	body := "" +
		"chr1\t101\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n" + // one read
		"chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n" + // no reads
		"chr2\t210\t.\tA\tC,G\t.\tPASS\t.\tGT\t1/2\n" // one read, two rows
	got, err := runVcfCount(t, body, "--min-depth", "1")
	if err != nil {
		t.Fatal(err)
	}
	want := "chrom\tpos\tref\talt\tref_count\talt_count\n" +
		"chr1\t101\tC\tT\t1\t0\n" +
		"chr2\t210\tA\tC\t1\t0\n" +
		"chr2\t210\tA\tG\t1\t0"
	if got != want {
		t.Errorf("--min-depth 1: got\n%s\nwant\n%s", got, want)
	}

	got, err = runVcfCount(t, body, "--min-depth", "2")
	if err != nil {
		t.Fatal(err)
	}
	if got != "chrom\tpos\tref\talt\tref_count\talt_count" {
		t.Errorf("--min-depth 2 should leave only the header:\n%s", got)
	}

	if _, err := runVcfCount(t, body, "--min-depth", "-1"); err == nil {
		t.Error("--min-depth -1 should be refused")
	}
}

func TestVcfCountMinDepthCountsOtherReads(t *testing.T) {
	// A site covered only by a third allele is covered: total is 1 though
	// ref and alt are both 0, so --min-depth 1 keeps it.
	v := vcVariantFor(t, 6, "C", "T")
	ref, alts, total := count(t, v, []*alignedRead{read("g", 0, 1, "10M", "ACGTAGGTAC", "")}, 13, true)
	if ref != 0 || alts[0] != 0 || total != 1 {
		t.Errorf("ref=%d alt=%d total=%d, want 0 0 1", ref, alts[0], total)
	}
}
