package vcfcmd

import (
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio/tabix"
)

// TestVcfOutputBgzipTbi covers -o NAME.gz writing BGZF (not plain gzip) and
// --tbi adding a queryable companion index.
func TestVcfOutputBgzipTbi(t *testing.T) {
	path := filepath.Join(t.TempDir(), "out.vcf.gz")
	runVcf(t, "vcf-concat", "-o", path, "--tbi", "testdata/concat_a.vcf", "testdata/concat_b.vcf")

	// BGZF blocks carry a "BC" extra subfield in the gzip header.
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.HasPrefix(data, []byte{0x1f, 0x8b}) || !bytes.Contains(data[:18], []byte("BC")) {
		t.Errorf("output is not BGZF-compressed: header = % x", data[:18])
	}

	if _, err := os.Stat(path + ".tbi"); err != nil {
		t.Fatalf("expected %s.tbi to be created: %v", path, err)
	}

	// The index must give random access to a record from either input file.
	r, err := tabix.NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	seq, err := r.Query("chr2", 499, 500)
	if err != nil {
		t.Fatal(err)
	}
	n := 0
	for rec, err := range seq {
		if err != nil {
			t.Fatal(err)
		}
		if !strings.HasPrefix(rec.Line, "chr2\t500\t") {
			t.Errorf("got %q", rec.Line)
		}
		n++
	}
	if n != 1 {
		t.Errorf("chr2:500 query returned %d rows, want 1", n)
	}
}

// TestVcfSplitTbi covers --tbi indexing every chunk vcf-split writes.
func TestVcfSplitTbi(t *testing.T) {
	base := filepath.Join(t.TempDir(), "out")
	runVcf(t, "vcf-split", "--out", base, "--num", "2", "--tbi", "testdata/sample.vcf")

	// sample.vcf holds 5 variants -> 3 chunks, each with its own index.
	for i := 1; i <= 3; i++ {
		chunk := base + "." + strconv.Itoa(i) + ".vcf.gz"
		if _, err := os.Stat(chunk + ".tbi"); err != nil {
			t.Errorf("expected %s.tbi: %v", chunk, err)
		}
	}
	if _, err := os.Stat(base + ".4.vcf.gz.tbi"); err == nil {
		t.Errorf("unexpected index for a 4th chunk")
	}

	// The first chunk's index must resolve one of its two variants.
	r, err := tabix.NewReader(base + ".1.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	seq, err := r.Query("chr1", 199, 200)
	if err != nil {
		t.Fatal(err)
	}
	n := 0
	for rec, err := range seq {
		if err != nil {
			t.Fatal(err)
		}
		if !strings.HasPrefix(rec.Line, "chr1\t200\t") {
			t.Errorf("got %q", rec.Line)
		}
		n++
	}
	if n != 1 {
		t.Errorf("chr1:200 query returned %d rows, want 1", n)
	}
}

// TestVcfTbiUnsorted covers the order guard: an unsorted output still gets
// written in full, but the index is refused instead of silently wrong.
func TestVcfTbiUnsorted(t *testing.T) {
	cases := []struct {
		name, input, want string
	}{
		{"backwards", "testdata/unsorted.vcf", "chr1:100, which follows chr1:300"},
		{"interleaved", "testdata/interleaved.vcf", "chr1:300, but chr1 already gave way to chr2"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "out.vcf.gz")
			err := runVcfErr(t, "vcf-filter", "-o", path, "--tbi", tc.input)
			if err == nil {
				t.Fatalf("expected an unsorted-output error")
			}
			if !strings.Contains(err.Error(), tc.want) {
				t.Errorf("error = %v, want it to name %s", err, tc.want)
			}

			// The VCF itself must be complete and readable; only the index is
			// skipped.
			if _, err := os.Stat(path + ".tbi"); err == nil {
				t.Errorf("index was written for unsorted output")
			}
			if n := countVcfRecords(t, path); n != 3 {
				t.Errorf("output holds %d records, want all 3", n)
			}
		})
	}
}

// TestVcfSplitTbiUnsorted covers the same guard per chunk in vcf-split.
func TestVcfSplitTbiUnsorted(t *testing.T) {
	base := filepath.Join(t.TempDir(), "out")
	err := runVcfErr(t, "vcf-split", "--out", base, "--num", "2", "--tbi", "testdata/unsorted.vcf")
	if err == nil {
		t.Fatalf("expected an unsorted-chunk error")
	}
	// Chunk 1 holds chr1:300 then chr1:100, so its index must be refused.
	if !strings.Contains(err.Error(), "chr1:100, which follows chr1:300") {
		t.Errorf("error = %v, want the out-of-order pair named", err)
	}
	if _, err := os.Stat(base + ".1.vcf.gz.tbi"); err == nil {
		t.Errorf("index was written for an unsorted chunk")
	}
	if n := countVcfRecords(t, base+".1.vcf.gz"); n != 2 {
		t.Errorf("chunk 1 holds %d records, want 2", n)
	}
}

// TestVcfTbiNeedsBgzipFile covers the --tbi preconditions: a real output file
// with a bgzip name.
func TestVcfTbiNeedsBgzipFile(t *testing.T) {
	if err := runVcfErr(t, "vcf-concat", "--tbi", "testdata/concat_a.vcf"); err == nil {
		t.Errorf("expected an error for --tbi writing to stdout")
	} else if !strings.Contains(err.Error(), "output file") {
		t.Errorf("error = %v, want output-file complaint", err)
	}

	plain := filepath.Join(t.TempDir(), "out.vcf")
	if err := runVcfErr(t, "vcf-concat", "-o", plain, "--tbi", "testdata/concat_a.vcf"); err == nil {
		t.Errorf("expected an error for --tbi with an uncompressed name")
	} else if !strings.Contains(err.Error(), ".gz") {
		t.Errorf("error = %v, want bgzip-name complaint", err)
	}
}
