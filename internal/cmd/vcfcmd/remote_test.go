package vcfcmd

import (
	"context"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// The rule that can regress silently. A URL cannot be stat'd, so classification
// has to test the scheme first -- but doing so must not change what happens to
// anything that is not a URL, or a mistyped locus starts reporting "no such
// file" instead of a locus error.
func TestTargetLocatorIsAFileNotASelector(t *testing.T) {
	for _, v := range []string{
		"https://host/panel.bed",
		"http://host/sites.txt",
		"s3://bucket/panel.vcf",
	} {
		if !isTargetFile(v) {
			t.Errorf("%q was classified as an inline selector, not a file", v)
		}
	}
}

func TestMistypedLocusStillGetsALocusError(t *testing.T) {
	// None of these contain "://", so they must take exactly the path they
	// always did: not a file, therefore parsed as a selector.
	for _, v := range []string{
		"chr1:100:A",
		"chr1:1O0-2000",
		"chr1:abc",
		"HLA-A*01:01:01:01",
		"chr1",
	} {
		if isTargetFile(v) {
			t.Errorf("%q was classified as a file", v)
		}
	}

	// And the error still names the locus rather than the filesystem.
	_, err := parseTargets(context.Background(), []string{"chr1:1O0-2000"})
	if err == nil {
		t.Fatal("a mistyped region parsed cleanly")
	}
	if strings.Contains(err.Error(), "no such file") {
		t.Errorf("a mistyped locus produced a file error: %v", err)
	}
}

// A contig name carrying colons still splits into four fields, and is still a
// contig, because its tail is not numeric. Pinned because the scheme check runs
// just before the code that decides this.
func TestColonBearingContigIsStillAContig(t *testing.T) {
	tset, err := parseTargets(context.Background(), []string{"HLA-A*01:01:01:01"})
	if err != nil {
		t.Fatalf("a colon-bearing contig name was rejected: %v", err)
	}
	if len(tset.Spans) != 1 || tset.Spans[0].Chrom != "HLA-A*01:01:01:01" {
		t.Errorf("got %+v, want one whole-contig span", tset.Spans)
	}
}

// Target files read over HTTP must be detected as the same format they would be
// locally. The sniff now reads through a SectionReader rather than an *os.File,
// and a difference there would flip BED-vs-site-list detection -- which shifts
// every coordinate by one rather than failing.
func TestTargetFileFormatIsTheSameOverHTTP(t *testing.T) {
	dir := t.TempDir()
	files := map[string]string{
		"panel.bed":  "chr1\t99\t100\nchr1\t199\t200\n",
		"sites.txt":  "chr1 100 A G\nchr1 200 C T\n",
		"panel.vcf":  "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t.\t.\t.\n",
		"tokens.txt": "chr1:100:A:G\n",
	}
	for name, body := range files {
		if err := os.WriteFile(filepath.Join(dir, name), []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	srv := httptest.NewServer(http.FileServer(http.Dir(dir)))
	defer srv.Close()

	ctx := context.Background()
	for name := range files {
		local, lerr := detectTargetFile(ctx, filepath.Join(dir, name))
		remote, rerr := detectTargetFile(ctx, srv.URL+"/"+name)
		if lerr != nil || rerr != nil {
			t.Fatalf("%s: local err %v, remote err %v", name, lerr, rerr)
		}
		if local != remote {
			t.Errorf("%s: detected as %q locally and %q over http", name, local, remote)
		}

		lt, err := parseTargets(ctx, []string{filepath.Join(dir, name)})
		if err != nil {
			t.Fatal(err)
		}
		rt, err := parseTargets(ctx, []string{srv.URL + "/" + name})
		if err != nil {
			t.Fatal(err)
		}
		if len(lt.Loci) != len(rt.Loci) || len(lt.Spans) != len(rt.Spans) {
			t.Errorf("%s: local gave %d loci/%d spans, remote %d/%d",
				name, len(lt.Loci), len(lt.Spans), len(rt.Loci), len(rt.Spans))
		}
		for i := range lt.Loci {
			if lt.Loci[i] != rt.Loci[i] {
				t.Errorf("%s locus %d: local %v, remote %v", name, i, lt.Loci[i], rt.Loci[i])
			}
		}
	}
}

func TestSampleFileFromLocator(t *testing.T) {
	dir := t.TempDir()
	if err := os.WriteFile(filepath.Join(dir, "subjects.txt"), []byte("S1\nS2\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	srv := httptest.NewServer(http.FileServer(http.Dir(dir)))
	defer srv.Close()

	set, err := parseSampleArgs(context.Background(), []string{srv.URL + "/subjects.txt"})
	if err != nil {
		t.Fatal(err)
	}
	if len(set.Names) != 2 || set.Names[0] != "S1" || set.Names[1] != "S2" {
		t.Errorf("got %v, want [S1 S2]", set.Names)
	}
}

// A remote input is read, a remote output is refused, and the refusal names the
// flag and says "local" rather than failing deep inside a writer.
func TestRemoteOutputIsRefused(t *testing.T) {
	for _, tc := range []struct {
		name string
		args []string
	}{
		{"toparquet --out", []string{"vcf-toparquet", "--out", "s3://bucket/cohort", "testdata/sample.vcf"}},
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
	err := runVcfErr(t, "vcf-varquery", "--variant", "chr1:100", "gs://bucket/cohort.vcf.gz")
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

// Reading a store over HTTP must give the same rows as reading it locally.
func TestStoreQueryOverHTTP(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "cohort")
	runVcf(t, "vcf-toparquet", "--out", base, "testdata/coverage.vcf")

	srv := httptest.NewServer(http.FileServer(http.Dir(dir)))
	defer srv.Close()

	local := runVcf(t, "vcf-varquery", "--variant", "1:100:A:G", base)
	remote := runVcf(t, "vcf-varquery", "--variant", "1:100:A:G", srv.URL+"/cohort")
	if local == "" {
		t.Fatal("the local query returned nothing; the fixture proves nothing")
	}
	if stripQueryProvenance(local) != stripQueryProvenance(remote) {
		t.Errorf("remote query differs from local\nlocal:\n%s\nremote:\n%s", local, remote)
	}
}

// stripQueryProvenance drops vcf-varquery's "## " banner, which echoes the input
// path and so differs by construction between a local and a remote read. The
// existing stripProvenance strips VCF header lines, which is a different set.
func stripQueryProvenance(s string) string {
	var keep []string
	for _, line := range strings.Split(s, "\n") {
		if strings.HasPrefix(line, "## ") {
			continue
		}
		keep = append(keep, line)
	}
	return strings.Join(keep, "\n")
}
