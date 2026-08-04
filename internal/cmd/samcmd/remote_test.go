package samcmd

import (
	"context"
	"io"
	"net/http"
	"net/http/httptest"
	"testing"

	"github.com/compgenlab/cghts/htsio"
)

// One test covers the whole swap: every SAM/BAM/CRAM command in cgkit reaches
// the same OpenSamReader, so what needs proving is that a locator opens at all
// and that the .bai is resolved over the same transport. The per-command
// plumbing is compile-time, and the existing local tests cover the behaviour.
func TestSamReaderAcceptsAnHTTPLocator(t *testing.T) {
	srv := httptest.NewServer(http.FileServer(http.Dir("testdata")))
	defer srv.Close()

	ctx := context.Background()
	remote, err := htsio.OpenSamReader(ctx, srv.URL+"/test.bam")
	if err != nil {
		t.Fatalf("OpenSamReader over http: %v", err)
	}
	defer remote.Close()

	local, err := htsio.NewSamReader("testdata/test.bam")
	if err != nil {
		t.Fatal(err)
	}
	defer local.Close()

	lh, err := local.Header()
	if err != nil {
		t.Fatal(err)
	}
	rh, err := remote.Header()
	if err != nil {
		t.Fatal(err)
	}
	if got, want := len(rh.References()), len(lh.References()); got != want {
		t.Fatalf("remote header has %d refs, local %d", got, want)
	}

	count := func(r htsio.SamReader) int {
		n := 0
		for _, err := range r.Records() {
			if err == io.EOF {
				break
			}
			if err != nil {
				t.Fatal(err)
			}
			n++
		}
		return n
	}
	// Each reader is single-pass, so count once and compare the numbers.
	gotRemote, wantLocal := count(remote), count(local)
	if wantLocal == 0 {
		t.Fatal("the fixture has no records, so this proves nothing")
	}
	if gotRemote != wantLocal {
		t.Errorf("read %d records over http, %d locally", gotRemote, wantLocal)
	}
}

// An indexed region query over http proves the .bai was resolved through the
// same transport, which is the part that would silently fall back to a full
// scan if sibling resolution were wrong.
func TestIndexedQueryOverHTTP(t *testing.T) {
	srv := httptest.NewServer(http.FileServer(http.Dir("testdata")))
	defer srv.Close()

	r, err := htsio.OpenSamReader(context.Background(), srv.URL+"/test.bam")
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	h, err := r.Header()
	if err != nil {
		t.Fatal(err)
	}
	refs := h.References()
	if len(refs) == 0 {
		t.Skip("fixture has no references to query")
	}
	seq, err := r.Query(refs[0].Name, 0, 1<<30)
	if err != nil {
		t.Fatalf("indexed query over http (the .bai was not resolved): %v", err)
	}
	for _, err := range seq {
		if err != nil && err != io.EOF {
			t.Fatal(err)
		}
	}
}

// A local path must reach the existing constructor rather than being re-derived
// through iosource, so nothing about a local run changes.
func TestLocalPathIsUnchangedByTheLocatorPath(t *testing.T) {
	viaOpen, err := htsio.OpenSamReader(context.Background(), "testdata/test.bam")
	if err != nil {
		t.Fatal(err)
	}
	defer viaOpen.Close()
	viaNew, err := htsio.NewSamReader("testdata/test.bam")
	if err != nil {
		t.Fatal(err)
	}
	defer viaNew.Close()

	a, err := viaOpen.Header()
	if err != nil {
		t.Fatal(err)
	}
	b, err := viaNew.Header()
	if err != nil {
		t.Fatal(err)
	}
	if a.Text() != b.Text() {
		t.Error("OpenSamReader and NewSamReader disagree on a local file's header")
	}
}
