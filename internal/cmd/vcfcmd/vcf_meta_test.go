package vcfcmd

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/varstore"
)

// The named flags are generated from varstore.ReservedMetaKeys, so every
// reserved key must reach the store under its own name. A drift here would put
// a value under the wrong key rather than fail, so it is checked per key.
func TestToParquetRecordsEveryReservedMetaKey(t *testing.T) {
	base := filepath.Join(t.TempDir(), "store")
	args := []string{"vcf-tovarstore", "--out", base}
	want := map[string]string{}
	for _, k := range varstore.ReservedMetaKeys {
		v := "value-for-" + k
		want[k] = v
		args = append(args, "--meta-"+k, v)
	}
	runVcf(t, append(args, "testdata/coverage.vcf")...)

	m, err := varstore.ReadManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	for k, v := range want {
		if m.Meta[k] != v {
			t.Errorf("manifest meta[%q] = %q, want %q", k, m.Meta[k], v)
		}
	}
	if len(m.Meta) != len(want) {
		t.Errorf("manifest meta has %d keys, want %d: %v", len(m.Meta), len(want), m.Meta)
	}
}

func TestToParquetGenericMetaFlag(t *testing.T) {
	base := filepath.Join(t.TempDir(), "store")
	runVcf(t, "vcf-tovarstore", "--out", base,
		"--meta-dataset", "1kg",
		"--meta", "cohort=phase3",
		"--meta", "batch=b7",
		// A value containing "=" keeps everything past the first separator.
		"--meta", "url=https://example.org/?a=1&b=2",
		"testdata/coverage.vcf")

	m, err := varstore.ReadManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	want := map[string]string{
		"dataset": "1kg",
		"cohort":  "phase3",
		"batch":   "b7",
		"url":     "https://example.org/?a=1&b=2",
	}
	for k, v := range want {
		if m.Meta[k] != v {
			t.Errorf("meta[%q] = %q, want %q", k, m.Meta[k], v)
		}
	}
}

// A key supplied twice is refused rather than resolved. Which of two
// conflicting claims about what a store holds gets recorded must not depend on
// flag order.
func TestToParquetRejectsDuplicateMetaKey(t *testing.T) {
	cases := []struct {
		name string
		args []string
	}{
		{"across spellings", []string{"--meta-dataset", "A", "--meta", "dataset=B"}},
		{"both generic", []string{"--meta", "cohort=A", "--meta", "cohort=B"}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			out := filepath.Join(t.TempDir(), "store")
			args := append([]string{"vcf-tovarstore", "--out", out}, c.args...)
			err := runVcfErr(t, append(args, "testdata/coverage.vcf")...)
			if err == nil {
				t.Fatal("expected an error for a duplicated metadata key")
			}
			if !strings.Contains(err.Error(), "given twice") {
				t.Errorf("error %q does not explain the duplication", err)
			}
		})
	}
}

func TestToParquetRejectsMalformedMeta(t *testing.T) {
	cases := []struct {
		name, arg, want string
	}{
		{"no separator", "nokey", "expected KEY=VALUE"},
		{"empty key", "=value", "empty key"},
		{"uppercase key", "Dataset=x", "invalid metadata key"},
		{"dotted key", "a.b=x", "invalid metadata key"},
		{"spaced key", "a b=x", "invalid metadata key"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			out := filepath.Join(t.TempDir(), "store")
			err := runVcfErr(t, "vcf-tovarstore", "--out", out,
				"--meta", c.arg, "testdata/coverage.vcf")
			if err == nil {
				t.Fatalf("expected an error for --meta %q", c.arg)
			}
			if !strings.Contains(err.Error(), c.want) {
				t.Errorf("error %q does not contain %q", err, c.want)
			}
			// Rejected before anything is created: the check runs ahead of
			// EnsureStoreDir so a typo leaves no directory to clean up.
			if _, statErr := os.Stat(out); statErr == nil {
				t.Errorf("a rejected conversion created %s", out)
			}
		})
	}
}

func TestVarSummaryReportsMetadata(t *testing.T) {
	base := filepath.Join(t.TempDir(), "store")
	runVcf(t, "vcf-tovarstore", "--out", base,
		"--meta-dataset", "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05",
		"--meta-reference", "GRCh38",
		"--meta", "cohort=phase3",
		"testdata/coverage.vcf")

	got := runVcf(t, "vcf-varsummary", base)
	if !strings.Contains(got, "metadata") {
		t.Fatalf("no metadata block in:\n%s", got)
	}
	for _, want := range []string{
		"dataset    20201028_CCDG_14151_B01_GRM_WGS_2020-08-05",
		"reference  GRCh38",
		"cohort     phase3",
	} {
		if !strings.Contains(got, want) {
			t.Errorf("summary missing %q:\n%s", want, got)
		}
	}
	// Reserved keys lead, in the library's order; unreserved keys follow.
	if strings.Index(got, "reference") > strings.Index(got, "cohort") {
		t.Errorf("reserved keys should precede unreserved ones:\n%s", got)
	}
}

// --format json emits the manifest verbatim, which is the promise that makes
// the open-ended map shape safe: anything recorded is reachable with jq even if
// the text report never learns to print it.
func TestVarSummaryJSONCarriesMetadata(t *testing.T) {
	base := filepath.Join(t.TempDir(), "store")
	runVcf(t, "vcf-tovarstore", "--out", base,
		"--meta-dataset", "1kg", "--meta", "cohort=phase3", "testdata/coverage.vcf")

	var doc struct {
		Meta map[string]string `json:"meta"`
	}
	if err := json.Unmarshal([]byte(runVcf(t, "vcf-varsummary", "--format", "json", base)), &doc); err != nil {
		t.Fatal(err)
	}
	if doc.Meta["dataset"] != "1kg" || doc.Meta["cohort"] != "phase3" {
		t.Errorf("json meta = %v", doc.Meta)
	}
}

// A store converted without metadata reports none, rather than an empty block.
// "Not stated" and "stated as nothing" are different claims.
func TestVarSummaryOmitsAbsentMetadata(t *testing.T) {
	base := convert(t, "testdata/coverage.vcf")

	if got := runVcf(t, "vcf-varsummary", base); strings.Contains(got, "metadata") {
		t.Errorf("metadata block present for a store with none:\n%s", got)
	}
	if got := runVcf(t, "vcf-varsummary", "--format", "json", base); strings.Contains(got, `"meta"`) {
		t.Errorf("json carries a meta key for a store with none:\n%s", got)
	}
}

// -v reports the two keys that change how the rows should be read, and only
// those; the rest is catalogue that belongs in vcf-varsummary.
func TestVarQueryVerboseReportsDatasetAndReference(t *testing.T) {
	base := filepath.Join(t.TempDir(), "store")
	runVcf(t, "vcf-tovarstore", "--out", base,
		"--meta-dataset", "1kg-phase3",
		"--meta-reference", "GRCh38",
		"--meta-caller", "GATK 4.2.6.1",
		"testdata/coverage.vcf")

	got := runVcf(t, "vcf-varquery", "-v", "--variant", "1:100:A:G", base)
	for _, want := range []string{"dataset     1kg-phase3", "reference   GRCh38"} {
		if !strings.Contains(got, want) {
			t.Errorf("verbose output missing %q:\n%s", want, got)
		}
	}
	if strings.Contains(got, "GATK") {
		t.Errorf("verbose output should not carry every metadata key:\n%s", got)
	}
}
