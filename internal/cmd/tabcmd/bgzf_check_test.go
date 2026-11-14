package tabcmd

import (
	"bytes"
	"compress/flate"
	"compress/gzip"
	"encoding/binary"
	"fmt"
	"hash/crc32"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio/bgzf"
	"github.com/compgenlab/cghts/htsio/tabix"
)

// bedLines is a sorted BED large enough to span several 64 KB BGZF blocks.
func bedLines() []string {
	var lines []string
	for _, chrom := range []string{"chr1", "chr2"} {
		for i := 0; i < 4000; i++ {
			lines = append(lines, fmt.Sprintf("%s\t%d\t%d\tfeature_%s_%d_with_some_padding_text", chrom, i*100, i*100+50, chrom, i))
		}
	}
	return lines
}

// writeBgzf writes lines the normal way, through cghts's BGZF writer.
func writeBgzf(t *testing.T, path string, lines []string) {
	t.Helper()
	f, err := os.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	w := bgzf.NewWriter(f)
	for _, l := range lines {
		if _, err := w.Write([]byte(l + "\n")); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}
}

// writeBlockPerLine writes each line as its own BGZF block -- valid BGZF, but
// the pathological layout bgzf-check exists to expose.
func writeBlockPerLine(t *testing.T, path string, lines []string, withEOF bool) {
	t.Helper()
	var out, comp bytes.Buffer
	fw, _ := flate.NewWriter(&comp, flate.DefaultCompression)
	for _, l := range lines {
		data := []byte(l + "\n")
		comp.Reset()
		fw.Reset(&comp)
		fw.Write(data)
		fw.Close()
		hdr := []byte{0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, 'B', 'C', 0x02, 0x00, 0, 0}
		binary.LittleEndian.PutUint16(hdr[16:], uint16(len(hdr)+comp.Len()+8-1))
		out.Write(hdr)
		out.Write(comp.Bytes())
		binary.Write(&out, binary.LittleEndian, crc32.ChecksumIEEE(data))
		binary.Write(&out, binary.LittleEndian, uint32(len(data)))
	}
	if withEOF {
		out.Write(bgzfEOFMarker)
	}
	if err := os.WriteFile(path, out.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
}

func indexBed(t *testing.T, path string) {
	t.Helper()
	if err := tabix.NewIndexWriter(tabix.NewWriterOpts().BED()).WriteIndex(path); err != nil {
		t.Fatalf("indexing %s: %v", path, err)
	}
}

func check(t *testing.T, path, index string, blocks bool) (string, error) {
	t.Helper()
	var out bytes.Buffer
	err := runBgzfCheck(&out, path, index, blocks)
	return out.String(), err
}

func TestBgzfCheckNormalFileAndIndex(t *testing.T) {
	path := filepath.Join(t.TempDir(), "ok.bed.gz")
	writeBgzf(t, path, bedLines())
	indexBed(t, path)

	got, err := check(t, path, "", false)
	if err != nil {
		t.Fatalf("a well-formed file and its own index should pass: %v\n%s", err, got)
	}
	for _, want := range []string{
		"EOF marker:\tpresent",
		"Small blocks (< 16.0 KB uncompressed, last block excluded):\t0 (0.0%)",
		"Index:\t" + path + ".tbi (TBI)",
		"Columns:\tseq=1 begin=2 end=3",
		"Zero-based:\ttrue",
		"Chunks not on a block boundary:\t0",
	} {
		if !strings.Contains(got, want) {
			t.Errorf("report lacks %q:\n%s", want, got)
		}
	}
	if strings.Contains(got, "[Warnings]") || strings.Contains(got, "[Problems]") {
		t.Errorf("unexpected warnings or problems:\n%s", got)
	}
}

func TestBgzfCheckSmallBlocksWarnButPass(t *testing.T) {
	path := filepath.Join(t.TempDir(), "small.bed.gz")
	writeBlockPerLine(t, path, bedLines(), true)
	indexBed(t, path)

	got, err := check(t, path, "", false)
	if err != nil {
		t.Fatalf("small blocks are valid BGZF and should only warn: %v\n%s", err, got)
	}
	if !strings.Contains(got, "Blocks:\t8000") {
		t.Errorf("want one block per line (8000):\n%s", got)
	}
	if !strings.Contains(got, "[Warnings]") || !strings.Contains(got, "bgzip fills blocks to ~64 KB") {
		t.Errorf("want a small-block warning:\n%s", got)
	}
	// Its own index matches it, offsets "just past" a block included.
	if !strings.Contains(got, "Chunks not on a block boundary:\t0") {
		t.Errorf("a small-block file's own index should match it:\n%s", got)
	}
}

func TestBgzfCheckMissingEOF(t *testing.T) {
	path := filepath.Join(t.TempDir(), "trunc.bed.gz")
	writeBlockPerLine(t, path, bedLines()[:10], false)
	got, err := check(t, path, "", false)
	if err == nil || !strings.Contains(got, "EOF marker:\tMISSING") {
		t.Errorf("err = %v, want a missing-EOF problem:\n%s", err, got)
	}
}

func TestBgzfCheckPlainGzipIsRefused(t *testing.T) {
	path := filepath.Join(t.TempDir(), "plain.bed.gz")
	var buf bytes.Buffer
	gz := gzip.NewWriter(&buf)
	gz.Write([]byte(strings.Join(bedLines(), "\n") + "\n"))
	gz.Close()
	if err := os.WriteFile(path, buf.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := check(t, path, "", false); err == nil || !strings.Contains(err.Error(), "not a valid BGZF") {
		t.Errorf("err = %v, want a not-BGZF error", err)
	}
}

// An index built for one layout of the data, applied to another layout of the
// same records, is the stale-index case: it parses fine and points nowhere.
func TestBgzfCheckMismatchedIndex(t *testing.T) {
	dir := t.TempDir()
	good := filepath.Join(dir, "good.bed.gz")
	other := filepath.Join(dir, "other.bed.gz")
	writeBgzf(t, good, bedLines())
	writeBlockPerLine(t, other, bedLines(), true)
	indexBed(t, other)

	got, err := check(t, good, other+".tbi", false)
	if err == nil {
		t.Fatalf("an index built for a different file should fail:\n%s", got)
	}
	if !strings.Contains(got, "not block boundaries of this file") {
		t.Errorf("want the block-boundary problem named:\n%s", got)
	}
}

func TestBgzfCheckBlocksDump(t *testing.T) {
	lines := bedLines()
	path := filepath.Join(t.TempDir(), "dump.bed.gz")
	writeBgzf(t, path, lines)
	got, err := check(t, path, "", true)
	if err != nil {
		t.Fatal(err)
	}
	rows := strings.Split(strings.TrimSpace(got), "\n")
	if rows[0] != "block\tcoffset\tcsize\tuoffset\tusize" || len(rows) < 3 {
		t.Fatalf("unexpected --blocks output:\n%s", got)
	}
	// Uncompressed sizes must add up to the content, and blocks must tile the
	// file: each starts where the previous ended.
	var usum, nextC int64
	for _, r := range rows[1:] {
		var i, coff, csize, uoff, usize int64
		fmt.Sscanf(r, "%d\t%d\t%d\t%d\t%d", &i, &coff, &csize, &uoff, &usize)
		if coff != nextC {
			t.Errorf("block %d starts at %d, want %d", i, coff, nextC)
		}
		nextC = coff + csize
		usum += usize
	}
	if want := int64(len(strings.Join(lines, "\n")) + 1); usum != want {
		t.Errorf("uncompressed total %d, want %d", usum, want)
	}
}

func TestBgzfCheckRefusesRemote(t *testing.T) {
	if _, err := check(t, "s3://bucket/file.gz", "", false); err == nil || !strings.Contains(err.Error(), "local files only") {
		t.Errorf("err = %v, want a local-only refusal", err)
	}
}

// A separately compressed header and body, concatenated, leaves an EOF marker
// (an empty block) in mid-file. cghts before v0.12.2 stopped there, so this
// check reported a one-block file and, with the index, blamed every offset past
// the join. Now the whole file is read, and the file's own index -- which the
// tabix writer could not build before either -- checks clean.
func TestBgzfCheckConcatenatedStreams(t *testing.T) {
	dir := t.TempDir()
	hdr := filepath.Join(dir, "hdr.gz")
	body := filepath.Join(dir, "body.gz")
	lines := bedLines()
	writeBgzf(t, hdr, []string{"#chrom\tstart\tend\tname"})
	writeBgzf(t, body, lines)
	a, _ := os.ReadFile(hdr)
	b, _ := os.ReadFile(body)
	path := filepath.Join(dir, "concat.bed.gz")
	if err := os.WriteFile(path, append(a, b...), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := tabix.NewIndexWriter(tabix.NewWriterOpts().BED().Meta('#')).WriteIndex(path); err != nil {
		t.Fatalf("indexing a concatenated file: %v", err)
	}

	got, err := check(t, path, "", false)
	if err != nil {
		t.Fatalf("a concatenated file with its own index should pass: %v\n%s", err, got)
	}
	want := int64(len("#chrom\tstart\tend\tname\n") + len(strings.Join(lines, "\n")) + 1)
	for _, w := range []string{
		"Uncompressed size:\t" + humanBytes(want),
		"Chunks not on a block boundary:\t0",
	} {
		if !strings.Contains(got, w) {
			t.Errorf("report lacks %q:\n%s", w, got)
		}
	}
	if strings.Contains(got, "reading stopped") {
		t.Errorf("must read through the mid-file marker:\n%s", got)
	}
}
