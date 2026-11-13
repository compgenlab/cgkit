package tabcmd

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/htsio/bgzf"
	"github.com/compgenlab/cghts/htsio/tabix"
	"github.com/compgenlab/cgkit/internal/locator"
	"github.com/spf13/cobra"
)

var bgzfCheckCmd = &cobra.Command{
	GroupID:     "tabcmd",
	Annotations: map[string]string{"since": "v0.6.11"},
	Use:         "bgzf-check <file.gz>",
	Short:       "Check a BGZF file's blocks and its tabix index",
	Long: `Check a BGZF-compressed file and the tabix index beside it, for troubleshooting
slow or memory-hungry indexed reads.

The file is read once, start to finish: every block is decompressed and its CRC
checked. The report gives:

  - whether the file is BGZF at all, and whether it ends with the EOF marker
  - the number of blocks, and min/median/mean/max of each block's compressed and
    uncompressed size. bgzip fills blocks to about 64 KB uncompressed; a file of
    much smaller blocks is valid but larger than it needs to be, and every
    indexed read of it costs more blocks.

When an index is found (<file>.tbi or <file>.csi, or --index), it is checked
against the blocks just read:

  - the indexed columns (seq/begin/end), comment character and coordinate base,
    so an index built with the wrong preset is visible
  - whether the index is older than the data file
  - every chunk's start and end offsets must land on a real block boundary of
    this file. An index built for a different version of the file usually
    fails here, and silently returns wrong or no records when queried.
  - the largest chunk, in compressed bytes: roughly the most a single query can
    be made to read

--blocks writes one row per block instead (block, coffset, csize, uoffset,
usize) for deeper digging.

Exits non-zero when it finds a problem: not BGZF, no EOF marker, or index
offsets that do not fit the file. Small blocks and an index older than the data
are reported as warnings only.`,
	Args: cobra.ExactArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		return runBgzfCheck(cmd.OutOrStdout(), args[0], bgzfCheckIndex, bgzfCheckBlocks)
	},
}

var (
	bgzfCheckIndex  string
	bgzfCheckBlocks bool
)

func init() {
	f := bgzfCheckCmd.Flags()
	f.StringVar(&bgzfCheckIndex, "index", "", "Index to check (default: <file>.tbi or <file>.csi if present)")
	f.BoolVar(&bgzfCheckBlocks, "blocks", false, "Write one row per block (block, coffset, csize, uoffset, usize) instead of the report")
}

// bgzfEOFMarker is the 28-byte empty block the SAM/BGZF spec requires at the
// end of every BGZF file. A file without it was most likely truncated.
var bgzfEOFMarker = []byte{
	0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
	0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
}

// smallBlockBytes is the uncompressed size below which a block counts as
// small. bgzip and cghts both fill blocks to ~64 KB; only a file's last block
// (and the block before an EOF-separated concatenation) is legitimately short.
const smallBlockBytes = 16 * 1024

// bgzfBlock is one data block: where it starts and how big it is, compressed
// and uncompressed. The EOF marker is not a data block.
type bgzfBlock struct {
	coff, csize int64
	uoff, usize int64
}

// bgzfScan is what one pass over a BGZF file finds.
type bgzfScan struct {
	size   int64 // file size on disk
	blocks []bgzfBlock
	total  int64 // uncompressed bytes
	hasEOF bool
}

// scanBgzf reads the whole file through cghts's IndexingReader, which
// decompresses and CRC-checks every block and records where each one starts.
func scanBgzf(path string) (*bgzfScan, error) {
	st, err := os.Stat(path)
	if err != nil {
		return nil, err
	}
	s := &bgzfScan{size: st.Size()}

	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	if st.Size() >= int64(len(bgzfEOFMarker)) {
		tail := make([]byte, len(bgzfEOFMarker))
		if _, err := f.ReadAt(tail, st.Size()-int64(len(tail))); err != nil {
			return nil, err
		}
		s.hasEOF = bytes.Equal(tail, bgzfEOFMarker)
	}

	ir := bgzf.NewIndexingReader(f)
	n, err := io.Copy(io.Discard, ir)
	if err != nil {
		return nil, fmt.Errorf("%s is not a valid BGZF file: %w", path, err)
	}
	s.total = n
	if st.Size() > 0 && n == 0 && !s.hasEOF {
		return nil, fmt.Errorf("%s is not a valid BGZF file: no blocks", path)
	}

	// Block 0 is implicit at (0,0); IndexingReader reports every block after it.
	starts := append([]bgzf.BlockOffset{{}}, ir.Blocks()...)
	dataEnd := s.size
	if s.hasEOF {
		dataEnd -= int64(len(bgzfEOFMarker))
	}
	if n == 0 {
		return s, nil // only the EOF marker: an empty file, validly compressed
	}
	for i, b := range starts {
		nextC, nextU := dataEnd, s.total
		if i+1 < len(starts) {
			nextC, nextU = starts[i+1].Compressed, starts[i+1].Uncompressed
		}
		s.blocks = append(s.blocks, bgzfBlock{
			coff: b.Compressed, csize: nextC - b.Compressed,
			uoff: b.Uncompressed, usize: nextU - b.Uncompressed,
		})
	}
	return s, nil
}

// blockAt returns the index of the block starting exactly at coff, or -1.
func (s *bgzfScan) blockAt(coff int64) int {
	i := sort.Search(len(s.blocks), func(i int) bool { return s.blocks[i].coff >= coff })
	if i < len(s.blocks) && s.blocks[i].coff == coff {
		return i
	}
	return -1
}

func runBgzfCheck(w io.Writer, path, indexPath string, dumpBlocks bool) error {
	if locator.IsRemote(path) {
		return fmt.Errorf("bgzf-check reads local files only: %s", path)
	}
	s, err := scanBgzf(path)
	if err != nil {
		return err
	}
	out := bufio.NewWriter(w)
	defer out.Flush()

	if dumpBlocks {
		fmt.Fprintln(out, "block\tcoffset\tcsize\tuoffset\tusize")
		for i, b := range s.blocks {
			fmt.Fprintf(out, "%d\t%d\t%d\t%d\t%d\n", i, b.coff, b.csize, b.uoff, b.usize)
		}
		return nil
	}

	var problems, warnings []string

	fmt.Fprintf(out, "File:\t%s\n", path)
	fmt.Fprintf(out, "Format:\tBGZF (all blocks decompressed and CRC-checked)\n")
	fmt.Fprintf(out, "Compressed size:\t%s\n", humanBytes(s.size))
	fmt.Fprintf(out, "Uncompressed size:\t%s\n", humanBytes(s.total))
	if s.total > 0 {
		fmt.Fprintf(out, "Compression ratio:\t%.2fx\n", float64(s.total)/float64(s.size))
	}
	if s.hasEOF {
		fmt.Fprintf(out, "EOF marker:\tpresent\n")
	} else {
		fmt.Fprintf(out, "EOF marker:\tMISSING\n")
		problems = append(problems, "no BGZF EOF marker: the file was probably truncated")
	}

	fmt.Fprintf(out, "\n[Blocks]\n")
	fmt.Fprintf(out, "Blocks:\t%d\n", len(s.blocks))
	if len(s.blocks) > 0 {
		cs := make([]int64, len(s.blocks))
		us := make([]int64, len(s.blocks))
		small := 0
		for i, b := range s.blocks {
			cs[i], us[i] = b.csize, b.usize
			// The last block is short by nature; it does not count as small.
			if b.usize < smallBlockBytes && i < len(s.blocks)-1 {
				small++
			}
		}
		fmt.Fprintf(out, "Compressed block size:\t%s\n", sizeStats(cs))
		fmt.Fprintf(out, "Uncompressed block size:\t%s\n", sizeStats(us))
		fmt.Fprintf(out, "Small blocks (< %s uncompressed, last block excluded):\t%d (%.1f%%)\n",
			humanBytes(smallBlockBytes), small, 100*float64(small)/float64(len(s.blocks)))
		if med := median(us); len(s.blocks) > 1 && med < smallBlockBytes {
			warnings = append(warnings, fmt.Sprintf(
				"median uncompressed block is %s; bgzip fills blocks to ~64 KB, so this file is larger than it needs to be and indexed reads touch more blocks",
				humanBytes(med)))
		}
	}

	idx, err := findIndex(path, indexPath)
	if err != nil {
		return err
	}
	if idx == "" {
		fmt.Fprintf(out, "\n[Index]\nIndex:\tnone found (looked for %s.tbi and %s.csi)\n", path, path)
	} else {
		p, wn, err := checkIndex(out, s, path, idx)
		if err != nil {
			return err
		}
		problems = append(problems, p...)
		warnings = append(warnings, wn...)
	}

	if len(warnings) > 0 {
		fmt.Fprintf(out, "\n[Warnings]\n")
		for _, m := range warnings {
			fmt.Fprintf(out, "- %s\n", m)
		}
	}
	if len(problems) > 0 {
		fmt.Fprintf(out, "\n[Problems]\n")
		for _, m := range problems {
			fmt.Fprintf(out, "- %s\n", m)
		}
		out.Flush()
		return fmt.Errorf("%d problem(s) found in %s", len(problems), path)
	}
	return nil
}

// findIndex returns the index to check: the one named, else <path>.tbi, else
// <path>.csi, else "".
func findIndex(path, named string) (string, error) {
	if named != "" {
		if _, err := os.Stat(named); err != nil {
			return "", fmt.Errorf("--index: %w", err)
		}
		return named, nil
	}
	for _, ext := range []string{".tbi", ".csi"} {
		if st, err := os.Stat(path + ext); err == nil && !st.IsDir() {
			return path + ext, nil
		}
	}
	return "", nil
}

// indexWindow is the region size checkIndex queries the index in, and
// maxIndexedPos the largest coordinate a tabix index can address (2^29).
const (
	indexWindow   = 1_000_000
	maxIndexedPos = 1 << 29
)

// tabixIndex is what checkIndex needs from either index format.
type tabixIndex interface {
	Query(refID int, start, end int) []tabix.Chunk
	RefNames() []string
}

type indexHeader struct {
	kind                                     string
	format, colSeq, colBeg, colEnd, meta, sk int32
	zeroBased                                bool
}

func loadIndex(path string) (tabixIndex, indexHeader, error) {
	if strings.HasSuffix(strings.ToLower(path), ".csi") {
		c, err := tabix.LoadCSI(path)
		if err != nil {
			return nil, indexHeader{}, fmt.Errorf("reading index %s: %w", path, err)
		}
		return c, indexHeader{"CSI", c.Format, c.ColSeq, c.ColBeg, c.ColEnd, c.Meta, c.Skip, c.ZeroBased}, nil
	}
	t, err := tabix.LoadTBI(path)
	if err != nil {
		return nil, indexHeader{}, fmt.Errorf("reading index %s: %w", path, err)
	}
	return t, indexHeader{"TBI", t.Format, t.ColSeq, t.ColBeg, t.ColEnd, t.Meta, t.Skip, t.ZeroBased}, nil
}

// checkIndex reports the index header and validates every chunk offset against
// the blocks of the data file.
func checkIndex(out io.Writer, s *bgzfScan, dataPath, idxPath string) (problems, warnings []string, err error) {
	idx, h, err := loadIndex(idxPath)
	if err != nil {
		return nil, nil, err
	}
	fmt.Fprintf(out, "\n[Index]\n")
	fmt.Fprintf(out, "Index:\t%s (%s)\n", idxPath, h.kind)

	if dst, e1 := os.Stat(dataPath); e1 == nil {
		if ist, e2 := os.Stat(idxPath); e2 == nil && ist.ModTime().Before(dst.ModTime()) {
			warnings = append(warnings, "the index is older than the data file; if the data was rewritten after indexing, the index no longer describes it")
		}
	}

	formats := map[int32]string{0: "generic", 1: "SAM", 2: "VCF"}
	fmt.Fprintf(out, "Format:\t%s\n", formats[h.format&0xffff])
	end := strconv.Itoa(int(h.colEnd))
	if h.colEnd == 0 {
		end = "0 (point features)"
	}
	fmt.Fprintf(out, "Columns:\tseq=%d begin=%d end=%s\n", h.colSeq, h.colBeg, end)
	meta := "none"
	if h.meta > 0 {
		meta = strconv.QuoteRune(rune(h.meta))
	}
	fmt.Fprintf(out, "Comment char:\t%s\n", meta)
	fmt.Fprintf(out, "Skip lines:\t%d\n", h.sk)
	fmt.Fprintf(out, "Zero-based:\t%v\n", h.zeroBased)
	names := idx.RefNames()
	fmt.Fprintf(out, "References:\t%d\n", len(names))

	// A virtual offset names a block start plus a byte within that block. The
	// byte may equal the block's size -- "just past its last byte", the same
	// position as the next block's start -- since index writers record a chunk
	// boundary as wherever the previous record ended. Offsets may also point at
	// the EOF marker or the end of the data.
	dataEnd := s.size
	if s.hasEOF {
		dataEnd -= int64(len(bgzfEOFMarker))
	}
	validOffset := func(vo bgzf.VirtualOffset) bool {
		c, u := vo.BlockOffset(), int64(vo.WithinBlock())
		if i := s.blockAt(c); i >= 0 {
			return u <= s.blocks[i].usize
		}
		return u == 0 && (c == dataEnd || c == s.size)
	}

	// The index is queried a window at a time rather than once per reference:
	// one whole-reference query returns a single merged chunk, which would check
	// two offsets per contig and say nothing about what a real query reads.
	// What one window query spans is also the useful size figure: the most a
	// region lookup can be made to read.
	type chunkKey struct{ b, e bgzf.VirtualOffset }
	seen := map[chunkKey]bool{}
	chunks, bad := 0, 0
	var examples []string
	var largest int64
	largestAt := ""
	for refID, name := range names {
		for start := 0; start < maxIndexedPos; start += indexWindow {
			var span int64
			for _, ch := range idx.Query(refID, start, start+indexWindow) {
				if okB, okE := validOffset(ch.Begin), validOffset(ch.End); !okB || !okE {
					if k := (chunkKey{ch.Begin, ch.End}); !seen[k] {
						seen[k] = true
						chunks++
						bad++
						if len(examples) < 5 {
							examples = append(examples, fmt.Sprintf("%s: chunk %d:%d-%d:%d",
								name, ch.Begin.BlockOffset(), ch.Begin.WithinBlock(), ch.End.BlockOffset(), ch.End.WithinBlock()))
						}
					}
					continue
				}
				if k := (chunkKey{ch.Begin, ch.End}); !seen[k] {
					seen[k] = true
					chunks++
				}
				span += ch.End.BlockOffset() - ch.Begin.BlockOffset()
			}
			if span > largest {
				largest = span
				largestAt = fmt.Sprintf("%s:%d-%d", name, start+1, start+indexWindow)
			}
		}
	}
	fmt.Fprintf(out, "Chunks checked:\t%d\n", chunks)
	fmt.Fprintf(out, "Chunks not on a block boundary:\t%d\n", bad)
	if largestAt != "" {
		fmt.Fprintf(out, "Largest read for a 1 Mb region:\t~%s compressed (%s)\n", humanBytes(largest), largestAt)
	}
	if bad > 0 {
		problems = append(problems, fmt.Sprintf(
			"%d of %d index chunks point at offsets that are not block boundaries of this file; the index was built for a different file (or a different version of this one). Examples: %s",
			bad, chunks, strings.Join(examples, "; ")))
	}
	return problems, warnings, nil
}

func sizeStats(v []int64) string {
	sorted := append([]int64(nil), v...)
	sort.Slice(sorted, func(i, j int) bool { return sorted[i] < sorted[j] })
	var sum int64
	for _, x := range sorted {
		sum += x
	}
	return fmt.Sprintf("min %s, median %s, mean %s, max %s",
		humanBytes(sorted[0]), humanBytes(median(sorted)),
		humanBytes(sum/int64(len(sorted))), humanBytes(sorted[len(sorted)-1]))
}

// median of v; v need not be sorted.
func median(v []int64) int64 {
	s := append([]int64(nil), v...)
	sort.Slice(s, func(i, j int) bool { return s[i] < s[j] })
	return s[len(s)/2]
}

func humanBytes(n int64) string {
	const k = 1024
	switch {
	case n < k:
		return fmt.Sprintf("%d B", n)
	case n < k*k:
		return fmt.Sprintf("%.1f KB", float64(n)/k)
	case n < k*k*k:
		return fmt.Sprintf("%.1f MB", float64(n)/(k*k))
	default:
		return fmt.Sprintf("%.2f GB", float64(n)/(k*k*k))
	}
}
