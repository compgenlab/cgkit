package tabcmd

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// writeGzip writes content to path as a gzip stream.
func writeGzip(t *testing.T, path, content string) {
	t.Helper()
	f, err := os.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	gw := gzip.NewWriter(f)
	if _, err := io.WriteString(gw, content); err != nil {
		t.Fatal(err)
	}
	if err := gw.Close(); err != nil {
		t.Fatal(err)
	}
}

// readBgzfLines returns the decompressed lines of a BGZF/gzip file (BGZF is a
// valid gzip stream, so the stdlib reader handles it).
func readBgzfLines(t *testing.T, path string) []string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	gr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	defer gr.Close()
	var lines []string
	sc := bufio.NewScanner(gr)
	for sc.Scan() {
		lines = append(lines, sc.Text())
	}
	if err := sc.Err(); err != nil {
		t.Fatal(err)
	}
	return lines
}

// records is already in the writer's output order (references in first-seen
// order, ascending within a reference), so these tests isolate the new input
// decompression and --skip pass-through rather than the writer's sort behavior.
var tabSortRecords = []string{
	"chr1\t90\t110\tgeneA",
	"chr1\t150\t160\tgeneB",
	"chr2\t400\t600\tgeneC",
}

func TestTabSortGzipInput(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.bed.gz")
	out := filepath.Join(dir, "out.bed.gz")

	// gzip-compressed input must be transparently decompressed (detected by
	// magic bytes); a clean round-trip proves the records were not read as raw
	// gzip bytes.
	writeGzip(t, in, strings.Join(tabSortRecords, "\n")+"\n")
	runTab(t, "tab-sort", "--preset", "bed", "--no-index", "-o", out, in)

	got := readBgzfLines(t, out)
	if strings.Join(got, "\n") != strings.Join(tabSortRecords, "\n") {
		t.Errorf("gzip input output:\n%q\nwant:\n%q", got, tabSortRecords)
	}
}

func TestTabSortPlainInput(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.bed")
	out := filepath.Join(dir, "out.bed.gz")

	// Plain (uncompressed) input must still work — the magic-byte check must not
	// disturb the non-gzip path.
	if err := os.WriteFile(in, []byte(strings.Join(tabSortRecords, "\n")+"\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	runTab(t, "tab-sort", "--preset", "bed", "--no-index", "-o", out, in)

	got := readBgzfLines(t, out)
	if strings.Join(got, "\n") != strings.Join(tabSortRecords, "\n") {
		t.Errorf("plain input output:\n%q\nwant:\n%q", got, tabSortRecords)
	}
}

func TestTabSortSkipHeader(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.bed")
	out := filepath.Join(dir, "out.bed.gz")

	// Two leading header lines must pass through verbatim, ahead of the data,
	// without being parsed as records (the second header has non-numeric
	// coordinate columns that would fail BED parsing if not skipped).
	header := []string{"# track header", "col1\tcol2\tcol3"}
	content := strings.Join(header, "\n") + "\n" + strings.Join(tabSortRecords, "\n") + "\n"
	if err := os.WriteFile(in, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	runTab(t, "tab-sort", "--preset", "bed", "--skip", "2", "--no-index", "-o", out, in)

	got := readBgzfLines(t, out)
	want := append(append([]string{}, header...), tabSortRecords...)
	if strings.Join(got, "\n") != strings.Join(want, "\n") {
		t.Errorf("skip-header output:\n%q\nwant:\n%q", got, want)
	}
}
