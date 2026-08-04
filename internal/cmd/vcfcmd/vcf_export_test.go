package vcfcmd

import (
	"strings"
	"testing"
)

// -1 is this file's sentinel for "every sample", and it was also what
// SampleIndex returned for a name it did not know. An unknown name therefore
// emitted one column in the header and one value per sample in each row: a
// table whose every column index is shifted, with nothing reported.
//
// "GT:9" covers the other half. SampleIndex resolves a numeric name
// positionally without a bounds check, so an out-of-range index used to pass
// straight through.
func TestExportRejectsAnUnknownFormatSample(t *testing.T) {
	for _, sel := range []string{"GT:NOSUCHSAMPLE", "GT:9"} {
		err := runVcfErr(t, "vcf-export", "--format", sel, "testdata/sample.vcf")
		if err == nil {
			t.Errorf("--format %s was accepted", sel)
			continue
		}
		if !strings.Contains(err.Error(), "no such sample") {
			t.Errorf("--format %s: unhelpful error %v", sel, err)
		}
	}
}

// dataColumns returns the column count of the header row and of the first data
// row, skipping the "##" provenance banner.
func dataColumns(out string) (header, first int) {
	for _, line := range strings.Split(strings.TrimSpace(out), "\n") {
		if line == "" || strings.HasPrefix(line, "##") {
			continue
		}
		n := len(strings.Split(line, "\t"))
		if header == 0 {
			header = n
			continue
		}
		return header, n
	}
	return header, 0
}

// The header and the data rows must agree on column count, which is precisely
// what the sentinel collision broke.
func TestExportColumnsLineUp(t *testing.T) {
	for _, sel := range []string{"GT", "GT:TUMOR"} {
		out := runVcf(t, "vcf-export", "--format", sel, "testdata/sample.vcf")
		h, f := dataColumns(out)
		if f == 0 {
			t.Errorf("--format %s produced no data rows:\n%s", sel, out)
			continue
		}
		if h != f {
			t.Errorf("--format %s: header has %d columns, a data row has %d\n%s", sel, h, f, out)
		}
	}
}

// A named sample selects only that sample's value; the unnamed form selects all
// of them. Both have to keep working, since the fix touches the branch that
// tells them apart.
func TestExportSampleSelection(t *testing.T) {
	oneH, _ := dataColumns(runVcf(t, "vcf-export", "--format", "GT:TUMOR", "testdata/sample.vcf"))
	allH, _ := dataColumns(runVcf(t, "vcf-export", "--format", "GT", "testdata/sample.vcf"))
	if allH <= oneH {
		t.Errorf("selecting one sample gave %d columns and selecting all gave %d; "+
			"all should be wider", oneH, allH)
	}
}

// selectByGlob filters --key against the header's declared FORMAT ids, so a key
// the header never declared yielded no column and said nothing -- and the "at
// least one field" guard has already passed by then, leaving a table of locus
// columns and nothing else. A warning, not an error: a glob matching nothing is
// reasonable across files with differing headers, and the keys that did match
// still export.
func TestSampleExportWarnsOnAnUndeclaredKey(t *testing.T) {
	warn := captureStderr(t, "vcf-sample-export", "--key", "NOSUCH", "testdata/sample.vcf")
	if !strings.Contains(warn, "NOSUCH") {
		t.Errorf("no warning naming the unmatched key:\n%s", warn)
	}

	// A key the header does declare must not warn.
	if w := captureStderr(t, "vcf-sample-export", "--key", "GT", "testdata/sample.vcf"); strings.Contains(w, "warning") {
		t.Errorf("warned about a declared key:\n%s", w)
	}

	// Nor must a glob that matched something.
	if w := captureStderr(t, "vcf-sample-export", "--key", "A*", "testdata/sample.vcf"); strings.Contains(w, "warning") {
		t.Errorf("warned about a glob that matched:\n%s", w)
	}
}
