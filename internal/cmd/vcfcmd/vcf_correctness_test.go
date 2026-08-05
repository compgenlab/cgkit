package vcfcmd

import (
	"strings"
	"testing"
)

// vcf-reorder used to warn about an unresolvable sample and carry on. Naming one
// sample wrong then produced a VCF one column short of the cohort that was asked
// for; naming them all wrong produced a header with no FORMAT and no sample
// columns over records that still had theirs -- not valid VCF at all, written
// with exit status 0.
func TestReorderRejectsUnknownSample(t *testing.T) {
	for _, args := range [][]string{
		{"vcf-reorder", "-s", "NORMAL,TYPO", "testdata/sample.vcf"},
		{"vcf-reorder", "-s", "NOPE,ALSONOPE", "testdata/sample.vcf"},
	} {
		err := runVcfErr(t, args...)
		if err == nil {
			t.Fatalf("%v: expected an error for an unknown sample", args)
		}
		if !strings.Contains(err.Error(), "no such sample") {
			t.Errorf("%v: error %q should name the problem", args, err)
		}
		// The roster is listed, because the usual cause is a typo and the fix is
		// visible only if the real names are in front of you.
		if !strings.Contains(err.Error(), "NORMAL") || !strings.Contains(err.Error(), "TUMOR") {
			t.Errorf("%v: error %q should list the samples the file has", args, err)
		}
	}
}

func TestReorderStillWorksWithValidSamples(t *testing.T) {
	got := runVcf(t, "vcf-reorder", "-s", "TUMOR,NORMAL", "testdata/sample.vcf")
	var header string
	for _, l := range strings.Split(got, "\n") {
		if strings.HasPrefix(l, "#CHROM") {
			header = l
			break
		}
	}
	if !strings.HasSuffix(header, "FORMAT\tTUMOR\tNORMAL") {
		t.Errorf("#CHROM line = %q, want the samples swapped", header)
	}
}

// header.SampleIndex maps a numeric name to n-1 with no upper bound, so
// "--sample 9" against a 2-sample file resolved to 8 and passed an `idx < 0`
// guard -- failing later, once per record, from inside the reader.
func TestOutOfRangeSampleIndexIsRejectedUpFront(t *testing.T) {
	cases := [][]string{
		{"vcf-tocount", "--sample", "9", "testdata/sample.vcf"},
		{"vcf-tobedpe", "--score", "DP:9", "testdata/sv.vcf"},
	}
	for _, args := range cases {
		err := runVcfErr(t, args...)
		if err == nil {
			t.Fatalf("%v: expected an error for an out-of-range sample index", args)
		}
		if !strings.Contains(err.Error(), "no such sample") {
			t.Errorf("%v: error %q should name the sample, not fail deep in the reader", args, err)
		}
		if strings.Contains(err.Error(), "out of range") {
			t.Errorf("%v: error %q is the late per-record failure this replaced", args, err)
		}
	}
}

// A command given no arguments must not write its help page into a redirected
// stdout and exit 0 -- "cgkit vcf-tobed > out.bed && next-step out.bed" then
// proceeds with a help page as its input.
func TestNoArgsIsAUsageErrorNotHelpOnStdout(t *testing.T) {
	for _, cmd := range []string{
		"vcf-tobed", "vcf-stats", "vcf-export", "vcf-samples", "vcf-check",
		"vcf-strip", "vcf-varquery", "vcf-varsummary", "vcf-toparquet",
	} {
		root, buf := vcfTestRoot(cmd)
		err := root.Execute()
		if err == nil {
			t.Errorf("%s with no arguments returned nil; a pipeline would treat that as success", cmd)
		}
		if strings.Contains(buf.String(), "Usage:") && err == nil {
			t.Errorf("%s printed usage but reported success", cmd)
		}
	}
}
