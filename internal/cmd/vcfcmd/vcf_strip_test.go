package vcfcmd

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// The help has always documented "@file"; the code used to probe with os.Stat
// instead, so the documented form did nothing and an undocumented bare filename
// worked. Both halves are pinned here.
func TestStripReadsNamesFromAtFile(t *testing.T) {
	dir := t.TempDir()
	list := filepath.Join(dir, "fields.txt")
	if err := os.WriteFile(list, []byte("DP\nAF\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	out := runVcf(t, "vcf-strip", "--info", "@"+list, "testdata/sample.vcf")
	for _, gone := range []string{"DP=", "AF="} {
		if strings.Contains(out, gone) {
			t.Errorf("@file did not strip %q:\n%s", gone, out)
		}
	}
	// Something not in the list must survive, or the test proves only that
	// everything was removed.
	if !strings.Contains(out, "##fileformat=VCF") {
		t.Errorf("output is not a VCF:\n%s", out)
	}
}

// A bare filename used to be read as a list. These values are short tokens like
// AC or DP, so that made behaviour depend on what happened to sit in the working
// directory. It is now an error naming the fix rather than a silent literal.
func TestStripBareFilenameIsAnError(t *testing.T) {
	dir := t.TempDir()
	list := filepath.Join(dir, "fields.txt")
	if err := os.WriteFile(list, []byte("DP\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	err := runVcfErr(t, "vcf-strip", "--info", list, "testdata/sample.vcf")
	if err == nil {
		t.Fatal("a bare filename was accepted as a field name")
	}
	if !strings.Contains(err.Error(), "@") {
		t.Errorf("the error does not say to use @: %v", err)
	}
}

// A value that is not a file stays a literal field name, which is the ordinary case.
func TestStripPlainNameIsLiteral(t *testing.T) {
	out := runVcf(t, "vcf-strip", "--info", "DP", "testdata/sample.vcf")
	if strings.Contains(out, "DP=") {
		t.Errorf("--info DP did not strip DP:\n%s", out)
	}
}
