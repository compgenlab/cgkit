package vcfcmd

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/varstore"
)

func TestVarsummaryOnStoreAndVcf(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	runVcf(t, "vcf-toparquet", "--out", base, "testdata/coverage.vcf")

	store := runVcf(t, "vcf-varsummary", base)
	vcf := runVcf(t, "vcf-varsummary", "testdata/coverage.vcf")

	// The sample roster is the one thing both backends know, and it must agree.
	storeSamples := runVcf(t, "vcf-varsummary", "--samples", base)
	vcfSamples := runVcf(t, "vcf-varsummary", "--samples", "testdata/coverage.vcf")
	if storeSamples != vcfSamples {
		t.Errorf("sample rosters differ\nstore:\n%s\nvcf:\n%s", storeSamples, vcfSamples)
	}
	if strings.TrimSpace(storeSamples) == "" {
		t.Fatal("no samples listed; the fixture proves nothing")
	}

	// Store-only facts must be reported as absent for a VCF, not fabricated.
	// Matched as a field at the start of a line: the VCF report mentions
	// --min-dp in prose, precisely to say it does not have one.
	hasField := func(out, field string) bool {
		for _, line := range strings.Split(out, "\n") {
			if strings.HasPrefix(line, field+" ") {
				return true
			}
		}
		return false
	}
	if !hasField(store, "min-dp") {
		t.Errorf("store summary omits the conversion --min-dp:\n%s", store)
	}
	if hasField(vcf, "min-dp") {
		t.Errorf("VCF summary reports a conversion --min-dp it cannot know:\n%s", vcf)
	}
	if hasField(vcf, "created") || hasField(vcf, "spans") {
		t.Errorf("VCF summary reports store-only provenance:\n%s", vcf)
	}
	if !strings.Contains(vcf, "provenance none") {
		t.Errorf("VCF summary does not say provenance is unavailable:\n%s", vcf)
	}
}

// The census counts rows that were written, which is what lets it contradict
// the input list and contig declarations stamped before the first record.
func TestVarsummaryCensusCountsRows(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	runVcf(t, "vcf-toparquet", "--out", base, "testdata/coverage.vcf")

	out := runVcf(t, "vcf-varsummary", "--counts", base)
	if !strings.Contains(out, "per chromosome") {
		t.Fatalf("--counts produced no census:\n%s", out)
	}

	m, err := varstore.ReadManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	if len(m.Chromosomes) == 0 {
		t.Fatal("manifest recorded no chromosomes")
	}
	for _, c := range m.Chromosomes {
		if !strings.Contains(out, c.Name) {
			t.Errorf("census omits %s:\n%s", c.Name, out)
		}
	}
}

// --format json emits the manifest verbatim. That is the promise that makes
// gzipping it acceptable: "| jq" is one pipe away rather than zero.
func TestVarsummaryJSONIsTheManifest(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	runVcf(t, "vcf-toparquet", "--out", base, "testdata/coverage.vcf")

	var fromCmd, fromDisk varstore.Manifest
	if err := json.Unmarshal([]byte(runVcf(t, "vcf-varsummary", "--format", "json", base)), &fromCmd); err != nil {
		t.Fatal(err)
	}
	m, err := varstore.ReadManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	fromDisk = *m

	if fromCmd.Counts != fromDisk.Counts {
		t.Errorf("counts differ: %+v vs %+v", fromCmd.Counts, fromDisk.Counts)
	}
	if len(fromCmd.Chromosomes) != len(fromDisk.Chromosomes) {
		t.Errorf("census lengths differ: %d vs %d", len(fromCmd.Chromosomes), len(fromDisk.Chromosomes))
	}
	if fromCmd.Params != fromDisk.Params {
		t.Errorf("params differ: %+v vs %+v", fromCmd.Params, fromDisk.Params)
	}
}

// The default report must be O(header) on both backends. A summary that
// silently reads a 200 GB VCF is worse than one that says what it will cost, and
// this is exactly the guarantee that erodes the first time somebody adds a field.
//
// The fixture is a VCF whose header is valid and whose records are not, so
// anything that reads past the header fails and anything that does not, does not.
func TestVarsummaryDefaultDoesNotScan(t *testing.T) {
	path := filepath.Join(t.TempDir(), "poisoned.vcf")
	body := strings.Join([]string{
		"##fileformat=VCFv4.2",
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
		"chr1\tNOT_A_POSITION\t.\tA\tT\t.\t.\t.\tGT\t0/1",
		"",
	}, "\n")
	if err := os.WriteFile(path, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}

	// The default report reads the header only, so the bad record is never seen.
	out := runVcf(t, "vcf-varsummary", path)
	if !strings.Contains(out, "samples") {
		t.Errorf("default report did not describe the file:\n%s", out)
	}

	// --counts reads records, so it must reach the bad one and fail.
	if err := runVcfErr(t, "vcf-varsummary", "--counts", path); err == nil {
		t.Error("--counts succeeded on a file whose records cannot be parsed, so it " +
			"did not actually read them")
	}
}

// A store missing its manifest cannot be queried, so vcf-varsummary has to be
// able to say why -- otherwise "no escape hatch" leaves a user with an
// unreadable store and no way to learn anything about it.
func TestVarsummaryExplainsAMissingManifest(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	runVcf(t, "vcf-toparquet", "--out", base, "testdata/coverage.vcf")
	if err := os.Remove(varstore.ManifestPath(base)); err != nil {
		t.Fatal(err)
	}

	err := runVcfErr(t, "vcf-varsummary", base)
	if err == nil {
		t.Fatal("a manifest-less store was summarized as though it were fine")
	}
	if !strings.Contains(err.Error(), varstore.ManifestFile) {
		t.Errorf("the error does not name what is missing: %v", err)
	}
}
