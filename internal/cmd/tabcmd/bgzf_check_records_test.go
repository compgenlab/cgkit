package tabcmd

import (
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio/bgzf"
	"github.com/compgenlab/cghts/htsio/tabix"
)

// bgzf-check reads each reference's record count from the index (its
// pseudo-bin) and checks the total against the records in the file.

func TestBgzfCheckRecordCounts(t *testing.T) {
	path := filepath.Join(t.TempDir(), "ok.bed.gz")
	writeBgzf(t, path, bedLines())
	indexBed(t, path)
	got, err := check(t, path, "", false)
	if err != nil {
		t.Fatalf("%v\n%s", err, got)
	}
	for _, want := range []string{
		"Records:\t8000 (from the index; matches the file)",
		"[Records per reference]\nchr1\t4000\nchr2\t4000\n",
	} {
		if !strings.Contains(got, want) {
			t.Errorf("report lacks %q:\n%s", want, got)
		}
	}
}

// Header lines are not records: comment lines (meta) and the first skip lines
// are excluded from the file's count, as tabix excludes them from the index's.
func TestBgzfCheckRecordCountsExcludeHeaderLines(t *testing.T) {
	dir := t.TempDir()
	lines := append([]string{"track name=x", "chrom\tstart\tend\tname", "#a comment"}, bedLines()...)
	lines = append(lines[:1000], append([]string{"#a comment between records"}, lines[1000:]...)...)
	path := filepath.Join(dir, "hdr.bed.gz")
	writeBgzf(t, path, lines)
	opts := tabix.NewWriterOpts().BED().Meta('#').Skip(2)
	if err := tabix.NewIndexWriter(opts).WriteIndex(path); err != nil {
		t.Fatal(err)
	}
	got, err := check(t, path, "", false)
	if err != nil || !strings.Contains(got, "Records:\t8000 (from the index; matches the file)") {
		t.Errorf("err = %v; want 8000 matching records:\n%s", err, got)
	}
}

// An index whose count disagrees with the file was built for another version
// of it. Here the records the index was built over are rewritten one short.
func TestBgzfCheckRecordCountMismatch(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "f.bed.gz")
	lines := bedLines()
	writeBlockPerLine(t, path, lines, true)
	indexBed(t, path)
	idx, err := os.ReadFile(path + ".tbi")
	if err != nil {
		t.Fatal(err)
	}
	// Drop the last record, keeping the index newer than the data.
	writeBlockPerLine(t, path, lines[:len(lines)-1], true)
	if err := os.WriteFile(path+".tbi", idx, 0o644); err != nil {
		t.Fatal(err)
	}
	got, err := check(t, path, "", false)
	if err == nil {
		t.Fatalf("a count mismatch must fail:\n%s", got)
	}
	for _, want := range []string{
		"Records:\t8000 in the index, 7999 in the file",
		"the index counts 8000 records but the file has 7999",
	} {
		if !strings.Contains(got, want) {
			t.Errorf("report lacks %q:\n%s", want, got)
		}
	}
}

// An index from before counts were written says nothing either way: reported,
// not a problem.
func TestBgzfCheckIndexWithoutCounts(t *testing.T) {
	path := filepath.Join(t.TempDir(), "old.bed.gz")
	writeBgzf(t, path, bedLines())
	indexBed(t, path)
	stripPseudoBins(t, path+".tbi")
	got, err := check(t, path, "", false)
	if err != nil {
		t.Fatalf("an index without counts is not a problem: %v\n%s", err, got)
	}
	if !strings.Contains(got, "Records:\tnot recorded in the index") || strings.Contains(got, "[Records per reference]") {
		t.Errorf("want counts reported as not recorded:\n%s", got)
	}
}

// htslib's index carries the same counts.
func TestBgzfCheckRecordCountsFromHtslib(t *testing.T) {
	bin, err := exec.LookPath("tabix")
	if err != nil {
		t.Skip("htslib tabix not on PATH")
	}
	path := filepath.Join(t.TempDir(), "hts.bed.gz")
	writeBgzf(t, path, bedLines())
	if out, err := exec.Command(bin, "-p", "bed", path).CombinedOutput(); err != nil {
		t.Fatalf("tabix: %v\n%s", err, out)
	}
	got, err := check(t, path, "", false)
	if err != nil || !strings.Contains(got, "Records:\t8000 (from the index; matches the file)") ||
		!strings.Contains(got, "chr1\t4000\nchr2\t4000\n") {
		t.Errorf("err = %v; want htslib's counts:\n%s", err, got)
	}
}

func TestLineCounter(t *testing.T) {
	data := "#h1\n\nskip me\nchr1\t1\r\n#c\nchr1\t2\nchr2\t3" // no final newline
	for _, size := range []int{1, 2, 3, 7, len(data)} {
		var c lineCounter
		for p := data; len(p) > 0; {
			n := min(size, len(p))
			c.Write([]byte(p[:n]))
			p = p[n:]
		}
		for _, tc := range []struct {
			meta, skip int32
			want       int64
		}{
			{0, 0, 6},   // every non-empty line
			{'#', 0, 4}, // comments out
			{'#', 3, 3}, // and the first three lines: "#h1", "", "skip me"
			{0, 3, 4},
		} {
			if got, ok := c.records(tc.meta, tc.skip); !ok || got != tc.want {
				t.Errorf("chunks of %d, meta %q skip %d: %d, %v; want %d", size, tc.meta, tc.skip, got, ok, tc.want)
			}
		}
	}
}

// stripPseudoBins rewrites a TBI as cghts wrote it before v0.12.9: without each
// reference's pseudo-bin (37450) or the trailing n_no_coor.
func stripPseudoBins(t *testing.T, path string) {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	zr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	d, err := io.ReadAll(zr)
	if err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	o := 0
	take := func(n int) []byte { b := d[o : o+n]; o += n; return b }
	i32 := func() int32 { return int32(binary.LittleEndian.Uint32(take(4))) }
	out.Write(take(4)) // magic
	nRef := i32()
	binary.Write(&out, binary.LittleEndian, nRef)
	out.Write(take(24)) // format .. skip
	lNm := i32()
	binary.Write(&out, binary.LittleEndian, lNm)
	out.Write(take(int(lNm)))
	for range nRef {
		nBin := i32()
		var bins bytes.Buffer
		kept := int32(0)
		for range nBin {
			bin := binary.LittleEndian.Uint32(take(4))
			nChunk := i32()
			chunks := take(16 * int(nChunk))
			if bin == 37450 {
				continue
			}
			kept++
			binary.Write(&bins, binary.LittleEndian, bin)
			binary.Write(&bins, binary.LittleEndian, nChunk)
			bins.Write(chunks)
		}
		binary.Write(&out, binary.LittleEndian, kept)
		out.Write(bins.Bytes())
		nIntv := i32()
		binary.Write(&out, binary.LittleEndian, nIntv)
		out.Write(take(8 * int(nIntv)))
	}
	if o != len(d)-8 {
		t.Fatalf("expected exactly n_no_coor after the references, %d bytes left", len(d)-o)
	}
	var z bytes.Buffer
	w := bgzf.NewWriter(&z)
	w.Write(out.Bytes())
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, z.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
}
