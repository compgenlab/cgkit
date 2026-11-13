package vcfcmd

import (
	"strings"
	"testing"
)

// A remote input is read, a remote output is refused, and the refusal names the
// flag and says "local" rather than failing deep inside a writer.
func TestRemoteOutputIsRefused(t *testing.T) {
	for _, tc := range []struct {
		name string
		args []string
	}{
		{"export -o", []string{"vcf-export", "-o", "https://host/out.tsv", "testdata/sample.vcf"}},
		{"tobed -o", []string{"vcf-tobed", "-o", "s3://bucket/out.bed", "testdata/sample.vcf"}},
	} {
		err := runVcfErr(t, tc.args...)
		if err == nil {
			t.Errorf("%s: a remote output was accepted", tc.name)
			continue
		}
		if !strings.Contains(err.Error(), "local") {
			t.Errorf("%s: error does not say output must be local: %v", tc.name, err)
		}
	}
}

// This test is what pins the s3 blank import to internal/locator. Move it to
// main.go and this still passes there while the shipped binary differs -- so the
// assertion is that the scheme is *registered*, from inside a package test.
func TestUnregisteredSchemeMessage(t *testing.T) {
	err := runVcfErr(t, "vcf-stats", "gs://bucket/cohort.vcf.gz")
	if err == nil {
		t.Fatal("a gs:// input was accepted")
	}
	msg := err.Error()
	if !strings.Contains(msg, "gs://") {
		t.Errorf("error does not name the scheme: %v", err)
	}
	if !strings.Contains(msg, "s3") {
		t.Errorf("error does not list s3 among the registered transports, so the "+
			"transport is not linked into this binary: %v", err)
	}
}
