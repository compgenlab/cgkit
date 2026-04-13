package ontcmd

import (
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/compgen-io/cgltk/htsio"
	"github.com/spf13/cobra"
)

var ontUmiClusterCmd = &cobra.Command{
	GroupID: "ontcmd",
	Use:     "ont-umi-cluster <input.bam>",
	Short:   "Collapse similar UMIs in a coordinate-sorted BAM file",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if umiClusterOutput == "" {
			return fmt.Errorf("--output is required")
		}
		if umiClusterRegion != "" && umiClusterWholeGenome {
			return fmt.Errorf("--region and --whole-genome are mutually exclusive")
		}

		inputFile := args[0]

		skipRefs := []string{}
		if umiClusterSkipRefs != "" {
			skipRefs = strings.Split(umiClusterSkipRefs, ",")
		}

		// Open counts writer if requested and emit the BED6+ header. The
		// header is a single '#'-prefixed line so downstream BED tools
		// skip it as a comment.
		var countsWriter io.Writer
		var closeCounts func() error
		if umiClusterCountsFilename != "" {
			var err error
			countsWriter, closeCounts, err = openWriter(umiClusterCountsFilename, true)
			if err != nil {
				return fmt.Errorf("opening umi-counts: %w", err)
			}
			defer func() {
				if err := closeCounts(); err != nil {
					fmt.Fprintf(os.Stderr, "error closing umi-counts: %v\n", err)
				}
			}()
		}

		if umiClusterWholeGenome {
			return umiClusterWholeGenomeMode(inputFile, skipRefs)
		}
		return umiClusterOverlapMode(inputFile, countsWriter, skipRefs)
	},
}

// bufferedRead holds a read in the overlap buffer along with precomputed
// coordinate fields and a union-find identifier.
type bufferedRead struct {
	id     int // global sequential ID for union-find
	rec    *htsio.SamRecord
	rname  string
	strand string // "+", "-", or "." for --no-strand
	start  int    // 0-based
	end    int    // 0-based, exclusive (start + CigarRefLen)
}

// unionFind implements a disjoint-set data structure with path compression
// and union by rank. The union method returns old/new roots so callers can
// merge associated maps.
type unionFind struct {
	parent []int
	rank   []int
}

func newUnionFind(capacity int) *unionFind {
	parent := make([]int, capacity)
	for i := range parent {
		parent[i] = i
	}
	return &unionFind{
		parent: parent,
		rank:   make([]int, capacity),
	}
}

// grow ensures the union-find has capacity for IDs [0, n).
func (uf *unionFind) grow(n int) {
	for len(uf.parent) < n {
		uf.parent = append(uf.parent, len(uf.parent))
		uf.rank = append(uf.rank, 0)
	}
}

func (uf *unionFind) find(x int) int {
	if uf.parent[x] != x {
		uf.parent[x] = uf.find(uf.parent[x])
	}
	return uf.parent[x]
}

// union merges the sets containing x and y. Returns the new root, old root,
// and whether a merge actually occurred (false if already in same set).
func (uf *unionFind) union(x, y int) (newRoot, oldRoot int, merged bool) {
	px, py := uf.find(x), uf.find(y)
	if px == py {
		return px, py, false
	}
	if uf.rank[px] < uf.rank[py] {
		px, py = py, px
	}
	uf.parent[py] = px
	if uf.rank[px] == uf.rank[py] {
		uf.rank[px]++
	}
	return px, py, true
}

// removeFromBin removes the bufferedRead with the given id from a bin
// slice using swap-with-last. Order within a bin doesn't matter, so this
// is O(bin_size) scan + O(1) removal. Returns the shortened slice.
func removeFromBin(bin []*bufferedRead, id int) []*bufferedRead {
	for i, r := range bin {
		if r.id == id {
			bin[i] = bin[len(bin)-1]
			bin[len(bin)-1] = nil
			return bin[:len(bin)-1]
		}
	}
	return bin
}


// umiClusterOverlapMode groups reads by 5' and/or 3' end proximity using a
// buffer + union-find approach, then clusters UMIs within each component.
//
// Parallelism model (rewritten 2026-04-11):
//
//   - When --region is set, cgltk processes that region only. The caller
//     is expected to be orchestrating per-region jobs externally (e.g.
//     one SLURM task per chromosome), so we do not also do per-chromosome
//     parallelism inside a single invocation. Skipped refs and unmapped
//     reads are *not* passed through in region mode — they are handled
//     by separate jobs.
//
//   - When --region is not set, we process each chromosome in the header
//     sequentially. The entire --threads budget is given to the per-group
//     UMI clustering step (clusterUMIs), which is the real
//     bottleneck. samtools sort is cheap and does not need many threads.
func umiClusterOverlapMode(inputFile string, countsWriter io.Writer, skipRefs []string) error {
	// Read header from the input file.
	hdrReader, err := htsio.NewSamReader(inputFile)
	if err != nil {
		return err
	}
	header, err := hdrReader.Header()
	if err != nil {
		return fmt.Errorf("failed to read header: %w", err)
	}
	hdrReader.Close()
	if header == nil {
		return fmt.Errorf("no header found in BAM file")
	}
	if err := validateCoordinateSorted(header); err != nil {
		return err
	}

	addUMIClusterPGLine(header)

	// Build skip set for fast lookup (only used when --region is not set).
	skipSet := make(map[string]bool)
	for _, r := range skipRefs {
		skipSet[r] = true
	}

	// Determine which regions to process. When --region is set we only
	// do that one region; otherwise iterate all references in header
	// order (minus the skip list).
	var regions []string
	if umiClusterRegion != "" {
		regions = []string{umiClusterRegion}
	} else {
		for _, ref := range header.References() {
			if skipSet[ref.Name] {
				continue
			}
			regions = append(regions, ref.Name)
		}
	}

	// Open writer — samtools sort handles merging output. samtools itself
	// uses very little CPU so we keep writer threads small and leave the
	// compute budget for clusterUMIs.
	wopts := htsio.SamWriterOptions(header).BAM().SortCoord().Threads(2)
	writer, err := htsio.NewSamWriter(umiClusterOutput, wopts)
	if err != nil {
		return err
	}

	var nextMI int64 = 1
	var totalReads int64
	var totalChanged int64

	// Open the input reader. In --region mode we query a specific
	// region; otherwise we read the entire file in one pass and let
	// processReads handle chromosome transitions, skip-ref pass-
	// through, and unmapped pass-through inline.
	var reader htsio.SamReader
	if umiClusterRegion != "" {
		ropts := htsio.NewSamReaderOpts().Region(umiClusterRegion)
		r, err := htsio.NewSamReader(inputFile, ropts)
		if err != nil {
			writer.Close()
			return err
		}
		reader = r
	} else {
		r, err := htsio.NewSamReader(inputFile)
		if err != nil {
			writer.Close()
			return err
		}
		reader = r
	}

	if err := processReads(reader, writer, skipSet,
		&nextMI, &totalReads, &totalChanged, countsWriter); err != nil {
		reader.Close()
		writer.Close()
		return err
	}
	reader.Close()

	fmt.Fprintf(os.Stderr, "Total reads: %d, UMIs corrected: %d\n",
		totalReads, totalChanged)
	return writer.Close()
}

// processReads runs the overlap-mode grouping + UMI clustering over a
// pre-opened SamReader. It handles chromosome transitions: when the
// RefName changes between consecutive records, the current buffer is
// flushed and the new chromosome starts with a clean slate.
//
// When skipSet is non-nil (full-file mode), records on skipped refs and
// unmapped records are written through unchanged without clustering.
// When skipSet is nil (--region mode), unmapped records are dropped
// (they're another job's concern).
func processReads(
	reader htsio.SamReader,
	writer *htsio.SamtoolsSamWriter,
	skipSet map[string]bool,
	nextMI *int64,
	totalReads *int64,
	totalChanged *int64,
	countsWriter io.Writer,
) error {

	// Per-region state. Reads are indexed by start and end position in
	// bin maps (bin size = overlap) so overlap queries are O(matches)
	// instead of O(buffer_size). The union-find and component maps are
	// the same as before; only the detection data structure changed.
	overlap := umiClusterOverlap
	if overlap <= 0 {
		overlap = 1 // safety: avoid division by zero in bin keys
	}

	startIndex := make(map[int][]*bufferedRead) // key = start / overlap
	endIndex := make(map[int][]*bufferedRead)   // key = end / overlap
	active := make(map[int]*bufferedRead)        // key = read id

	uf := newUnionFind(1024)
	activeCount := make(map[int]int)
	componentReads := make(map[int][]*bufferedRead)
	globalID := 0

	// loStartBin / loEndBin are watermarks for the lowest bin that might
	// still contain reads. Ejection walks forward from the watermark
	// instead of scanning all bins.
	loStartBin := 0
	loEndBin := 0

	// lastEjectStart caches the most recent curStart passed to
	// ejectExpired so we can skip redundant calls when consecutive
	// reads share the same start position (common in high-depth regions).
	lastEjectStart := -1

	// -----------------------------------------------------------------
	// Worker pool for parallel cluster processing.
	//
	// The main goroutine handles detection (streaming reads, indexed
	// bins, union-find). When a cluster is complete (all its reads have
	// been ejected from the active set), it's sent to workCh. A pool
	// of worker goroutines pull clusters from workCh and run the
	// expensive part — clusterUMIs (all-pairs Levenshtein) — in
	// parallel. Each worker also writes its finished records to the
	// SamWriter (which is internally thread-safe via a buffered
	// channel) and to countsWriter (serialized via ioMu).
	//
	// Shared state that's mutated by workers:
	//   - nextMI: unique MI counter (atomic.Int64 for lock-free increment)
	//   - totalReads/totalChanged: summary counters (atomic.Int64)
	//   - countsWriter: serialized via ioMu
	//   - stderr log lines: serialized via ioMu
	//   - writeErr: first worker error (atomic.Value)
	//
	// The expensive clusterUMIs call itself can use multiple threads
	// internally (the all-pairs loop is already parallelized via
	// sync.WaitGroup inside clusterUMIs). When only one large cluster
	// is in flight, it gets the full CPU; when many small clusters are
	// in flight, the workers saturate cores from the pool level. Go's
	// runtime scheduler handles the multiplexing.
	// -----------------------------------------------------------------

	var ioMu sync.Mutex     // serializes countsWriter, stderr, and counts-related I/O
	var writeErr atomic.Value // first error from a worker (nil or error)

	var atomicNextMI atomic.Int64
	atomicNextMI.Store(*nextMI)
	var atomicTotalReads atomic.Int64
	var atomicTotalChanged atomic.Int64

	workCh := make(chan []*bufferedRead, umiClusterThreads*2)
	var workerWg sync.WaitGroup

	for i := 0; i < max(umiClusterThreads, 1); i++ {
		workerWg.Add(1)
		go func() {
			defer workerWg.Done()
			for reads := range workCh {
				if writeErr.Load() != nil {
					continue // drain channel but skip work
				}
				if len(reads) == 0 {
					continue
				}

				// Phase 1 (parallel-safe): extract coordinates and
				// cluster UMIs. This is the expensive O(unique²) step
				// that benefits from running on multiple workers.
				minStart := reads[0].start
				maxEnd := reads[0].end
				rname := reads[0].rname
				strand := reads[0].strand
				recs := make([]*htsio.SamRecord, len(reads))
				for i, br := range reads {
					recs[i] = br.rec
					if br.start < minStart {
						minStart = br.start
					}
					if br.end > maxEnd {
						maxEnd = br.end
					}
				}

				umiCounts := make(map[string]int)
				for _, rec := range recs {
					if umi := getUMI(rec); umi != "" {
						umiCounts[umi]++
					}
				}
				representative := make(map[string]string)
				results := clusterUMIs(umiCounts, representative, umiClusterThreads)

				// Build per-UMI coordinate bounding boxes BEFORE
				// updateRecordUMI rewrites the tags — the map must be
				// keyed by the original UMI string (which is what
				// umiClusterResult.umi contains), not the normalized
				// representative that updateRecordUMI writes back.
				// We use the bufferedRead's pre-computed start/end to
				// avoid re-parsing CIGARs.
				var coordsByUMI map[string]umiCoords
				if countsWriter != nil {
					coordsByUMI = make(map[string]umiCoords, len(reads))
					for _, br := range reads {
						umi := getUMI(br.rec)
						if umi == "" {
							continue
						}
						c, ok := coordsByUMI[umi]
						if !ok {
							coordsByUMI[umi] = umiCoords{minStart: br.start, maxEnd: br.end}
						} else {
							if br.start < c.minStart {
								c.minStart = br.start
							}
							if br.end > c.maxEnd {
								c.maxEnd = br.end
							}
							coordsByUMI[umi] = c
						}
					}
				}

				// Update UMI tags on records (each worker owns its own
				// records, so no lock needed for the tag mutations).
				changed := 0
				for _, rec := range recs {
					if updateRecordUMI(rec, representative) {
						changed++
					}
				}

				// Phase 2 (serialized via ioMu): assign MI values, build
				// counts lines, print the log line. These touch shared
				// counters and writers so they need the lock, but they're
				// fast relative to the clustering above.
				representativeCount, maxClustSize := clusterStats(results)
				strandLabel := ""
				if !umiClusterNoStrand {
					strandLabel = "(" + strand + ") "
				}

				repToMI := make(map[string]string)

				ioMu.Lock()
				if umiClusterMI || countsWriter != nil {
					for _, r := range results {
						if _, ok := repToMI[r.representative]; !ok {
							mi := atomicNextMI.Add(1) - 1
							repToMI[r.representative] = fmt.Sprintf("mi_%09d", mi)
						}
					}
				}

				fmt.Fprintf(os.Stderr, "%s:%d-%d: %s%d reads, %d unique UMIs -> %d representative (max cluster: %d)\n",
					rname, minStart, maxEnd, strandLabel,
					len(recs), len(results), representativeCount, maxClustSize)
				ioMu.Unlock()

				// Set MI tags on records now that repToMI is populated.
				if umiClusterMI {
					for _, rec := range recs {
						origUMI := getUMI(rec)
						if origUMI != "" {
							if rep, ok := representative[origUMI]; ok {
								if mi, ok := repToMI[rep]; ok {
									rec.Tags["MI"] = htsio.SamTag{Type: 'Z', Value: mi}
								}
							}
						}
					}
				}

				// Phase 3: write records and counts. SamWriter is
				// thread-safe (internal channel), countsWriter is not.
				var firstErr error
				for _, rec := range recs {
					if err := writer.Write(rec); err != nil {
						firstErr = fmt.Errorf("writing record to output BAM: %w", err)
						break
					}
				}

				if firstErr == nil && countsWriter != nil && len(results) > 0 {
					cl := buildUMIClusterCounts(rname, strand, results, repToMI, coordsByUMI)
					ioMu.Lock()
					for _, line := range cl {
						if _, err := fmt.Fprintln(countsWriter, line); err != nil {
							firstErr = fmt.Errorf("writing umi-counts line: %w", err)
							break
						}
					}
					ioMu.Unlock()
				}

				atomicTotalReads.Add(int64(len(recs)))
				atomicTotalChanged.Add(int64(changed))

				if firstErr != nil {
					writeErr.CompareAndSwap(nil, firstErr)
				}
			}
		}()
	}

	// submitComponent sends a completed union-find component to the
	// worker pool for clustering + writing.
	submitComponent := func(root int) {
		reads := componentReads[root]
		delete(componentReads, root)
		delete(activeCount, root)
		if writeErr.Load() != nil {
			return
		}
		workCh <- reads
	}

	mergeComponents := func(newRoot, oldRoot int) {
		activeCount[newRoot] += activeCount[oldRoot]
		delete(activeCount, oldRoot)
		if reads, ok := componentReads[oldRoot]; ok {
			componentReads[newRoot] = append(componentReads[newRoot], reads...)
			delete(componentReads, oldRoot)
		}
	}

	// ejectExpired walks the low-numbered bins and ejects reads that
	// have fallen behind the current position.
	//
	// Default mode: eject reads with start < curStart - overlap.
	//   → walk startIndex bins from loStartBin upward.
	// match-one-end mode: eject reads with end < curStart - overlap.
	//   → walk endIndex bins from loEndBin upward.
	//
	// Complexity: O(number of reads ejected) amortized, because each
	// read is ejected exactly once and the watermark only moves forward.
	ejectExpired := func(curStart int) {
		threshold := curStart - overlap
		threshBin := threshold / overlap

		if umiClusterMatchOneEnd {
			// Eject by end position: reads whose end is too far
			// behind can no longer match any future read's 3'-end,
			// and future reads' 5'-starts have already passed them.
			for bin := loEndBin; bin < threshBin; bin++ {
				for _, b := range endIndex[bin] {
					// Remove from startIndex + active only (endIndex
					// bin will be deleted wholesale below).
					sBin := b.start / overlap
					startIndex[sBin] = removeFromBin(startIndex[sBin], b.id)
					if len(startIndex[sBin]) == 0 {
						delete(startIndex, sBin)
					}
					delete(active, b.id)
					// Union-find bookkeeping.
					root := uf.find(b.id)
					componentReads[root] = append(componentReads[root], b)
					activeCount[root]--
					if activeCount[root] == 0 {
						submitComponent(root)
					}
				}
				delete(endIndex, bin)
			}
			// Partial bin at the threshold boundary: check individual reads.
			if reads, ok := endIndex[threshBin]; ok {
				kept := 0
				for _, b := range reads {
					if b.end < threshold {
						sBin := b.start / overlap
						startIndex[sBin] = removeFromBin(startIndex[sBin], b.id)
						if len(startIndex[sBin]) == 0 {
							delete(startIndex, sBin)
						}
						delete(active, b.id)
						root := uf.find(b.id)
						componentReads[root] = append(componentReads[root], b)
						activeCount[root]--
						if activeCount[root] == 0 {
							submitComponent(root)
						}
					} else {
						reads[kept] = b
						kept++
					}
				}
				for i := kept; i < len(reads); i++ {
					reads[i] = nil
				}
				if kept == 0 {
					delete(endIndex, threshBin)
				} else {
					endIndex[threshBin] = reads[:kept]
				}
			}
			loEndBin = threshBin
		} else {
			// Default mode: eject by start position. Since BAM is
			// coordinate-sorted, all reads in low start-bins have
			// starts too far behind curStart to match any future
			// read's 5'-end — and the buffer only contains reads
			// whose 5' check passes.
			for bin := loStartBin; bin < threshBin; bin++ {
				for _, b := range startIndex[bin] {
					eBin := b.end / overlap
					endIndex[eBin] = removeFromBin(endIndex[eBin], b.id)
					if len(endIndex[eBin]) == 0 {
						delete(endIndex, eBin)
					}
					delete(active, b.id)
					root := uf.find(b.id)
					componentReads[root] = append(componentReads[root], b)
					activeCount[root]--
					if activeCount[root] == 0 {
						submitComponent(root)
					}
				}
				delete(startIndex, bin)
			}
			// Partial bin at threshold boundary.
			if reads, ok := startIndex[threshBin]; ok {
				kept := 0
				for _, b := range reads {
					if b.start < threshold {
						eBin := b.end / overlap
						endIndex[eBin] = removeFromBin(endIndex[eBin], b.id)
						if len(endIndex[eBin]) == 0 {
							delete(endIndex, eBin)
						}
						delete(active, b.id)
						root := uf.find(b.id)
						componentReads[root] = append(componentReads[root], b)
						activeCount[root]--
						if activeCount[root] == 0 {
							submitComponent(root)
						}
					} else {
						reads[kept] = b
						kept++
					}
				}
				for i := kept; i < len(reads); i++ {
					reads[i] = nil
				}
				if kept == 0 {
					delete(startIndex, threshBin)
				} else {
					startIndex[threshBin] = reads[:kept]
				}
			}
			loStartBin = threshBin
		}
	}

	// ejectAll drains every remaining active read into their components
	// and submits completed components to the worker pool.
	ejectAll := func() {
		for _, b := range active {
			if writeErr.Load() != nil {
				return
			}
			root := uf.find(b.id)
			componentReads[root] = append(componentReads[root], b)
			activeCount[root]--
			if activeCount[root] == 0 {
				submitComponent(root)
			}
		}
		// Wipe all indices.
		for k := range startIndex {
			delete(startIndex, k)
		}
		for k := range endIndex {
			delete(endIndex, k)
		}
		for k := range active {
			delete(active, k)
		}
	}

	readStrand := func(rec *htsio.SamRecord) string {
		if umiClusterNoStrand {
			return "."
		}
		if rec.IsReverse() {
			return "-"
		}
		return "+"
	}

	// currentChrom tracks the reference name so we can detect chromosome
	// transitions. When the chromosome changes, we flush the entire
	// active buffer (different chromosomes can't share coordinates).
	currentChrom := ""

	for {
		if e := writeErr.Load(); e != nil {
			close(workCh)
			workerWg.Wait()
			return e.(error)
		}
		rec, err := reader.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			close(workCh)
			workerWg.Wait()
			return err
		}

		// Unmapped reads: in full-file mode (skipSet non-nil) pass
		// through unchanged; in --region mode drop them.
		if rec.IsUnmapped() || rec.Cigar == "*" {
			if skipSet != nil {
				if err := writer.Write(rec); err != nil {
					close(workCh)
					workerWg.Wait()
					return fmt.Errorf("writing unmapped record: %w", err)
				}
			}
			continue
		}

		// Skipped refs: pass through unchanged in full-file mode.
		if skipSet != nil && skipSet[rec.RefName] {
			if err := writer.Write(rec); err != nil {
				close(workCh)
				workerWg.Wait()
				return fmt.Errorf("writing skipped-ref record: %w", err)
			}
			continue
		}

		// Chromosome transition: flush the active buffer. Reads on
		// different chromosomes occupy different coordinate spaces and
		// cannot overlap, so everything in the buffer belongs to the
		// old chromosome and must be submitted before we start the
		// new one. Also reset the bin watermarks.
		if rec.RefName != currentChrom {
			if currentChrom != "" {
				ejectAll()
				if e := writeErr.Load(); e != nil {
					close(workCh)
					workerWg.Wait()
					return e.(error)
				}
			}
			currentChrom = rec.RefName
			loStartBin = 0
			loEndBin = 0
			lastEjectStart = -1
			fmt.Fprintf(os.Stderr, "Processing %s...\n", currentChrom)
		}

		readStart := rec.Pos - 1
		readEnd := readStart + htsio.CigarRefLen(rec.Cigar)
		strand := readStrand(rec)

		// Skip redundant ejection when curStart hasn't changed. Many
		// consecutive reads share the same start in high-depth regions;
		// the watermark already advanced on the first call, so repeated
		// calls with the same threshold just re-scan the boundary bin
		// for nothing.
		if readStart != lastEjectStart {
			ejectExpired(readStart)
			lastEjectStart = readStart
		}

		br := &bufferedRead{
			id:     globalID,
			rec:    rec,
			rname:  rec.RefName,
			strand: strand,
			start:  readStart,
			end:    readEnd,
		}
		uf.grow(globalID + 1)
		globalID++

		activeCount[br.id] = 1

		// Query the bin indices for overlapping reads. Track myRoot so
		// we can skip candidates that are already in our component —
		// once a mega-component has formed, nearly every bin entry
		// shares the same root, and skipping them avoids millions of
		// redundant strand/coordinate checks.
		myRoot := br.id

		endLo := (br.end - overlap) / overlap
		endHi := (br.end + overlap) / overlap

		if umiClusterMatchOneEnd {
			// 3'-end matches.
			for bin := endLo; bin <= endHi; bin++ {
				for _, b := range endIndex[bin] {
					if uf.find(b.id) == myRoot {
						continue // already in our component
					}
					if b.strand != br.strand {
						continue
					}
					diff := br.end - b.end
					if diff < 0 {
						diff = -diff
					}
					if diff <= overlap {
						newRoot, oldRoot, merged := uf.union(br.id, b.id)
						if merged {
							mergeComponents(newRoot, oldRoot)
							myRoot = newRoot
						}
					}
				}
			}
			// 5'-start matches.
			startLo := (br.start - overlap) / overlap
			startHi := br.start / overlap
			for bin := startLo; bin <= startHi; bin++ {
				for _, b := range startIndex[bin] {
					if uf.find(b.id) == myRoot {
						continue
					}
					if b.strand != br.strand {
						continue
					}
					if br.start-b.start <= overlap {
						newRoot, oldRoot, merged := uf.union(br.id, b.id)
						if merged {
							mergeComponents(newRoot, oldRoot)
							myRoot = newRoot
						}
					}
				}
			}
		} else {
			// Default mode (match-both-ends): all active reads
			// already pass the 5' check (start-based ejection
			// enforces it). Only the 3'-end proximity matters.
			for bin := endLo; bin <= endHi; bin++ {
				for _, b := range endIndex[bin] {
					if uf.find(b.id) == myRoot {
						continue
					}
					if b.strand != br.strand {
						continue
					}
					diff := br.end - b.end
					if diff < 0 {
						diff = -diff
					}
					if diff <= overlap {
						newRoot, oldRoot, merged := uf.union(br.id, b.id)
						if merged {
							mergeComponents(newRoot, oldRoot)
							myRoot = newRoot
						}
					}
				}
			}
		}

		// Insert into indices.
		sBin := br.start / overlap
		startIndex[sBin] = append(startIndex[sBin], br)
		eBin := br.end / overlap
		endIndex[eBin] = append(endIndex[eBin], br)
		active[br.id] = br
	}

	ejectAll()

	// Signal workers that no more clusters are coming, then wait for
	// all in-flight clusters to finish processing + writing.
	close(workCh)
	workerWg.Wait()

	// Propagate atomic counters back to the caller's plain int64 pointers.
	// This is safe because workers are done by this point.
	*nextMI = atomicNextMI.Load()
	*totalReads += atomicTotalReads.Load()
	*totalChanged += atomicTotalChanged.Load()

	if e := writeErr.Load(); e != nil {
		return e.(error)
	}

	if currentChrom != "" {
		fmt.Fprintf(os.Stderr, "Finished %s\n", currentChrom)
	}
	return nil
}

// umiClusterWholeGenomeMode uses two passes over the entire file:
// pass 1 collects all UMI counts, pass 2 rewrites with representative UMIs.
func umiClusterWholeGenomeMode(inputFile string, skipRefs []string) error {
	// Pass 1: collect all UMIs
	reader, err := htsio.NewSamReader(inputFile)
	if err != nil {
		return err
	}
	umiCounts := make(map[string]int)

	header, err := reader.Header()
	if err != nil {
		return fmt.Errorf("failed to read header: %w", err)
	}
	if header == nil {
		return fmt.Errorf("no header found in BAM file")
	}

	// Build a set for O(1) skip-ref lookups.
	skipSet := make(map[string]bool, len(skipRefs))
	for _, r := range skipRefs {
		skipSet[r] = true
	}

	// Pass 1: collect UMI counts from mapped reads on non-skipped refs.
	// Unmapped reads have no position and shouldn't influence whole-
	// genome clustering; skipped refs are passed through unchanged in
	// pass 2.
	for {
		rec, err := reader.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if rec.IsUnmapped() || rec.Cigar == "*" {
			continue
		}
		if skipSet[rec.RefName] {
			continue
		}
		if umi := getUMI(rec); umi != "" {
			umiCounts[umi]++
		}
	}
	reader.Close()

	// Cluster
	globalRepresentative := make(map[string]string)
	results := clusterUMIs(umiCounts, globalRepresentative, umiClusterThreads)
	representativeCount, maxClustSize := clusterStats(results)
	fmt.Fprintf(os.Stderr, "whole-genome: %d unique UMIs -> %d representative (max cluster: %d)\n", len(umiCounts), representativeCount, maxClustSize)

	// Pass 2: rewrite BAM
	addUMIClusterPGLine(header)

	reader2, err := htsio.NewSamReader(inputFile)
	if err != nil {
		return err
	}

	// Keep writer threads consistent with overlap mode: always pass 2 so
	// samtools sort can compress in parallel when it wants to. The cost
	// at --threads 1 is negligible.
	opts := htsio.SamWriterOptions(header).BAM().SortCoord().Threads(2)

	writer, err := htsio.NewSamWriter(umiClusterOutput, opts)
	if err != nil {
		reader2.Close()
		return err
	}

	// Pass 2: rewrite BAM. Unmapped reads and skipped-ref reads are
	// passed through unchanged; everything else has its UMI tag
	// rewritten to the cluster representative.
	//
	// Every error return below this point must close both reader2 and
	// writer to avoid leaking an orphan samtools sort child process.
	changed := 0
	total := 0
	for {
		rec, err := reader2.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			reader2.Close()
			writer.Close()
			return err
		}
		total++

		if rec.IsUnmapped() || rec.Cigar == "*" {
			if err := writer.Write(rec); err != nil {
				reader2.Close()
				writer.Close()
				return fmt.Errorf("writing unmapped record: %w", err)
			}
			continue
		}
		if skipSet[rec.RefName] {
			if err := writer.Write(rec); err != nil {
				reader2.Close()
				writer.Close()
				return fmt.Errorf("writing skipped-ref record: %w", err)
			}
			continue
		}

		if updateRecordUMI(rec, globalRepresentative) {
			changed++
		}

		if err := writer.Write(rec); err != nil {
			reader2.Close()
			writer.Close()
			return fmt.Errorf("writing record: %w", err)
		}
	}
	reader2.Close()

	fmt.Fprintf(os.Stderr, "Total reads: %d, UMIs corrected: %d\n", total, changed)
	return writer.Close()
}

// umiCoords holds the bounding box of all reads carrying a specific UMI.
type umiCoords struct {
	minStart int
	maxEnd   int
}

// buildUMIClusterCounts returns tab-delimited BED6+ counts lines for one
// read-overlap-group. Coordinates are per-cluster (the bounding box of
// all reads whose UMIs belong to that cluster), not per-component.
func buildUMIClusterCounts(
	rname, strand string,
	results []umiClusterResult,
	repToMI map[string]string,
	coordsByUMI map[string]umiCoords,
) []string {
	// Group results by representative UMI and compute per-cluster
	// bounding boxes by merging the coordinates of each member UMI.
	type clusterInfo struct {
		mi             string
		representative string
		numReads       int
		umis           []string
		maxEditDist    int
		minStart       int
		maxEnd         int
	}
	clusters := make(map[string]*clusterInfo)
	for _, r := range results {
		ci, ok := clusters[r.representative]
		if !ok {
			ci = &clusterInfo{
				mi:             repToMI[r.representative],
				representative: r.representative,
				minStart:       1<<63 - 1, // MaxInt
				maxEnd:         0,
			}
			clusters[r.representative] = ci
		}
		ci.numReads += r.count
		ci.umis = append(ci.umis, r.umi)
		if r.maxIntraClustDist > ci.maxEditDist {
			ci.maxEditDist = r.maxIntraClustDist
		}
		// Merge this UMI's read coordinates into the cluster's bounding box.
		if coords, ok := coordsByUMI[r.umi]; ok {
			if coords.minStart < ci.minStart {
				ci.minStart = coords.minStart
			}
			if coords.maxEnd > ci.maxEnd {
				ci.maxEnd = coords.maxEnd
			}
		}
	}

	// Output format is BED6+ (standard BED6 columns followed by extras):
	//   chrom, start, end, name, score, strand,
	//   representative, numUMIs, maxEditDist, umis
	// Coordinates are per-cluster (not per-component), so each cluster's
	// BED interval reflects where its reads actually map.
	var lines []string
	for _, ci := range clusters {
		lines = append(lines, fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%s",
			rname, ci.minStart, ci.maxEnd, ci.mi, ci.numReads, strand,
			ci.representative, len(ci.umis), ci.maxEditDist,
			strings.Join(ci.umis, ",")))
	}
	return lines
}


func getUMI(rec *htsio.SamRecord) string {
	tag, ok := rec.Tags[umiClusterTag]
	if !ok {
		return ""
	}
	return tag.Value
}

func validateCoordinateSorted(header *htsio.SamHeader) error {
	for _, line := range header.Lines {
		if strings.HasPrefix(line, "@HD\t") {
			if strings.Contains(line, "SO:coordinate") {
				return nil
			}
			return fmt.Errorf("BAM file is not coordinate-sorted (expected SO:coordinate in @HD header)")
		}
	}
	return fmt.Errorf("BAM file has no @HD header line; cannot verify sort order")
}

// detectSeparator returns the UMI separator: "/", "-", or "TT".
func detectSeparator(umi string) string {
	if strings.Contains(umi, "/") {
		return "/"
	}
	if strings.Contains(umi, "-") {
		return "-"
	}
	return "TT"
}

// normalizeUMISeparator converts all UMI separator formats to use "/" as the
// canonical separator. "/" is safe to pass directly to the MSA aligner since
// it is not the alignment gap character ("-").
func normalizeUMISeparator(umi string) string {
	switch detectSeparator(umi) {
	case "-":
		return strings.ReplaceAll(umi, "-", "/")
	case "TT":
		return strings.ReplaceAll(umi, "TT", "/")
	default:
		return umi
	}
}

// updateRecordUMI rewrites the UMI tag to its representative value if one exists.
// Returns true if the tag was changed.
func updateRecordUMI(rec *htsio.SamRecord, representative map[string]string) bool {
	tag, hasTag := rec.Tags[umiClusterTag]
	if !hasTag {
		return false
	}
	umi := tag.Value
	cons, ok := representative[umi]
	if !ok || cons == umi {
		return false
	}
	rec.Tags[umiClusterOrigTag] = htsio.SamTag{Type: 'Z', Value: umi}
	rec.Tags[umiClusterTag] = htsio.SamTag{Type: 'Z', Value: cons}
	return true
}

// addUMIClusterPGLine appends the ont-umi-cluster @PG header line.
func addUMIClusterPGLine(header *htsio.SamHeader) {
	header.AddPGLine("ont-umi-cluster", "cgltk", fmt.Sprintf("DS:UMI collapsing; representative UMI written to %s, original preserved in %s, molecule ID in MI", umiClusterTag, umiClusterOrigTag))
}

// levBuf holds reusable DP buffers for Levenshtein computation.
type levBuf struct {
	prev []int
	curr []int
}

// levDist computes the Levenshtein edit distance between a and b, optionally
// bounded by maxDist. A negative maxDist means "no bound" — compute the full
// distance.
//
// When maxDist >= 0, the function uses Ukkonen's cutoff: after filling each
// DP row it checks the row minimum, and if that minimum already exceeds
// maxDist it returns maxDist+1 immediately. The correctness argument is the
// monotonicity invariant "min(row i+1) >= min(row i)" — once the row
// minimum climbs above the bound, no subsequent row can drop back below
// it, so the final distance is provably > maxDist.
//
// This bound makes a massive difference for UMI clustering where the
// edit-distance threshold is small (e.g. 3) and the vast majority of
// pairs are far apart: instead of a full 20x20 = 400-cell DP, we exit
// after 3–4 rows on most pairs.
//
// When the distance actually exceeds the bound, the returned value is
// `maxDist + 1` (a sentinel meaning "greater than maxDist") — callers
// that just need the ≤ maxDist check can compare `dist <= maxDist`, and
// callers that want an exact-but-capped number can clamp the return to
// maxDist.
func levDist(a, b string, buf *levBuf, maxDist int) int {
	m, n := len(a), len(b)
	if m == 0 {
		return n
	}
	if n == 0 {
		return m
	}
	if cap(buf.prev) <= n {
		buf.prev = make([]int, n+1)
		buf.curr = make([]int, n+1)
	} else {
		buf.prev = buf.prev[:n+1]
		buf.curr = buf.curr[:n+1]
	}
	for j := 0; j <= n; j++ {
		buf.prev[j] = j
	}
	for i := 1; i <= m; i++ {
		buf.curr[0] = i
		rowMin := i
		for j := 1; j <= n; j++ {
			if a[i-1] == b[j-1] {
				buf.curr[j] = buf.prev[j-1]
			} else {
				buf.curr[j] = 1 + min(buf.prev[j], buf.curr[j-1], buf.prev[j-1])
			}
			if buf.curr[j] < rowMin {
				rowMin = buf.curr[j]
			}
		}
		// Ukkonen cutoff: if the best value in this row already exceeds
		// maxDist, no cell in any subsequent row can drop back below it.
		if maxDist >= 0 && rowMin > maxDist {
			return maxDist + 1
		}
		buf.prev, buf.curr = buf.curr, buf.prev
	}
	return buf.prev[n]
}

// levDistHP computes an HP-aware Levenshtein edit distance. It uses the
// same bounded DP as levDist, but insertions and deletions within a
// homopolymer run cost 0 instead of 1. Specifically:
//
//   - Deleting a[i] costs 0 if a[i] == a[i-1] (shrinking an HP run in a)
//   - Inserting b[j] costs 0 if b[j] == b[j-1] (shrinking an HP run in b)
//   - Substitutions always cost 1 (or 0 for a match)
//
// This correctly handles ONT's most common error mode (HP-length
// variation) without the distortion that full HP compression causes.
// For example:
//
//	AAGA vs AAAA: substitution G→A = distance 1 (correct, not inflated)
//	AAACG vs AACG: HP indel = distance 0 (HP error collapsed)
//	AAGA vs AAGGA: HP indel in G run = distance 0 (HP error collapsed)
func levDistHP(a, b string, buf *levBuf, maxDist int) int {
	m, n := len(a), len(b)
	if m == 0 {
		return n
	}
	if n == 0 {
		return m
	}
	if cap(buf.prev) <= n {
		buf.prev = make([]int, n+1)
		buf.curr = make([]int, n+1)
	} else {
		buf.prev = buf.prev[:n+1]
		buf.curr = buf.curr[:n+1]
	}
	// Initialize row 0: inserting b[0..j-1]. HP indels in b cost 0.
	buf.prev[0] = 0
	for j := 1; j <= n; j++ {
		cost := 1
		if j >= 2 && b[j-1] == b[j-2] {
			cost = 0 // HP indel in b
		}
		buf.prev[j] = buf.prev[j-1] + cost
	}
	for i := 1; i <= m; i++ {
		// Cost of deleting a[i-1] (row header).
		delCostA := 1
		if i >= 2 && a[i-1] == a[i-2] {
			delCostA = 0 // HP indel in a
		}
		buf.curr[0] = buf.prev[0] + delCostA
		rowMin := buf.curr[0]
		for j := 1; j <= n; j++ {
			if a[i-1] == b[j-1] {
				// Match: no cost.
				buf.curr[j] = buf.prev[j-1]
			} else {
				// Deletion cost: removing a[i-1]. Free if it extends
				// an HP run in a (a[i-1] == a[i-2]).
				del := buf.prev[j] + delCostA
				// Insertion cost: inserting b[j-1]. Free if it extends
				// an HP run in b (b[j-1] == b[j-2]).
				insCostB := 1
				if j >= 2 && b[j-1] == b[j-2] {
					insCostB = 0
				}
				ins := buf.curr[j-1] + insCostB
				// Substitution: always costs 1.
				sub := buf.prev[j-1] + 1
				buf.curr[j] = min(del, ins, sub)
			}
			if buf.curr[j] < rowMin {
				rowMin = buf.curr[j]
			}
		}
		if maxDist >= 0 && rowMin > maxDist {
			return maxDist + 1
		}
		buf.prev, buf.curr = buf.curr, buf.prev
	}
	return buf.prev[n]
}

// collisionProb returns the probability that two independent random UMIs
// of length L (over a 4-letter alphabet) are within Levenshtein edit
// distance d of each other, using the substitution-only approximation:
//
//	P(dist ≤ d) ≈ Σ_{k=0}^{d} C(L,k) × 3^k / 4^L
//
// This is an upper bound (indels expand the neighborhood slightly) but
// a good approximation for short UMIs where most errors are substitutions
// or short HP indels.
func collisionProb(L, d int) float64 {
	total := 0.0
	choose := 1.0 // C(L, k) computed iteratively
	pow3 := 1.0   // 3^k
	pow4L := math.Pow(4.0, float64(L))
	for k := 0; k <= d; k++ {
		total += choose * pow3 / pow4L
		// Update for next k: C(L, k+1) = C(L, k) * (L-k) / (k+1)
		choose *= float64(L-k) / float64(k+1)
		pow3 *= 3.0
	}
	return total
}

// computeRepresentativeUMI picks the representative UMI for a cluster.
// The most common UMI (by read count) is chosen. Ties are broken first by
// longer normalized length and then lexicographically so the choice is
// deterministic across runs.
func computeRepresentativeUMI(members []umiCount) string {
	if len(members) == 0 {
		return ""
	}
	if len(members) == 1 {
		return normalizeUMISeparator(members[0].umi)
	}

	best := 0
	bestNorm := normalizeUMISeparator(members[0].umi)
	for i := 1; i < len(members); i++ {
		norm := normalizeUMISeparator(members[i].umi)
		if members[i].count > members[best].count ||
			(members[i].count == members[best].count && len(norm) > len(bestNorm)) ||
			(members[i].count == members[best].count && len(norm) == len(bestNorm) && norm < bestNorm) {
			best = i
			bestNorm = norm
		}
	}
	return bestNorm
}

// umiClusterResult holds the clustering result for one UMI.
type umiClusterResult struct {
	umi               string
	representative    string
	count             int
	maxIntraClustDist int // maximum pairwise edit distance between any two members of this cluster
}

type umiCount struct {
	umi   string
	count int
}

// umiEdge represents a pair of UMIs within the edit distance threshold.
type umiEdge struct{ i, j, dist int }

// clusterUMIs clusters a set of UMI strings using all-pairs Levenshtein
// edit distance followed by single-linkage union-find. Pairs within
// umiClusterEditThreshold become edges; connected components become
// clusters. For each cluster, a representative UMI is chosen (the most
// common one; see computeRepresentativeUMI) and recorded in
// globalRepresentative (orig UMI -> representative UMI).
//
// numThreads controls the parallelism of the all-pairs loop. Pass 1 for a
// serial run; pass umiClusterThreads (or another positive value) to split
// the work across workers.
func clusterUMIs(umiCounts map[string]int, globalRepresentative map[string]string, numThreads int) []umiClusterResult {
	if len(umiCounts) <= 1 {
		// Single UMI or empty, nothing to cluster. Normalize the
		// representative so single-member clusters produce the same
		// separator encoding as multi-member clusters (the multi-UMI
		// path below runs every representative through
		// computeRepresentativeUMI → normalizeUMISeparator).
		var results []umiClusterResult
		for umi, count := range umiCounts {
			norm := normalizeUMISeparator(umi)
			globalRepresentative[umi] = norm
			results = append(results, umiClusterResult{umi: umi, representative: norm, count: count})
		}
		return results
	}

	if numThreads <= 0 {
		numThreads = 1
	}

	// Sort for stable ordering; count does not affect cluster membership.
	umis := make([]umiCount, 0, len(umiCounts))
	for umi, count := range umiCounts {
		umis = append(umis, umiCount{umi, count})
	}
	sort.Slice(umis, func(i, j int) bool {
		if umis[i].count != umis[j].count {
			return umis[i].count > umis[j].count
		}
		return umis[i].umi < umis[j].umi
	})

	n := len(umis)

	// Pre-normalize all UMI strings once.
	normalized := make([]string, n)
	for i, u := range umis {
		normalized[i] = normalizeUMISeparator(u.umi)
	}

	// Select the distance function. When --hp-dist is set, use an
	// HP-aware Levenshtein where insertions/deletions within a
	// homopolymer run cost 0 instead of 1. This correctly handles
	// ONT's most common error mode (HP-length variation) without the
	// distortion that full HP compression causes on non-HP mutations.
	distFn := levDist
	if umiClusterHPDist {
		distFn = levDistHP
	}

	// Determine the effective edit distance threshold. When
	// --adaptive-threshold is set, reduce the threshold for large
	// components where random collisions at higher distances become
	// likely. The collision probability for two independent UMIs of
	// length L at distance ≤ d is:
	//   P(d) ≈ Σ_{k=0}^{d} C(L,k) × 3^k / 4^L
	// The expected number of false pairs at distance ≤ d in a
	// component of N unique UMIs is:
	//   E = N×(N-1)/2 × P(d)
	// We reduce the threshold until E < maxFalsePairs (default 1).
	edgeThreshold := umiClusterEditThreshold
	if umiClusterAdaptiveThreshold && n > 1 {
		// Auto-detect UMI length from the first normalized string,
		// excluding separator characters.
		umiLen := 0
		for _, c := range normalized[0] {
			if c != '/' {
				umiLen++
			}
		}
		if umiLen > 0 {
			nPairs := float64(n) * float64(n-1) / 2.0
			for d := edgeThreshold; d >= 1; d-- {
				p := collisionProb(umiLen, d)
				expectedFalse := nPairs * p
				if expectedFalse > umiClusterMaxFalsePairs {
					edgeThreshold = d - 1
				}
			}
			if edgeThreshold < umiClusterEditThreshold {
				fmt.Fprintf(os.Stderr, "  adaptive threshold: %d -> %d (N=%d, expected false at %d: %.1f)\n",
					umiClusterEditThreshold, edgeThreshold, n,
					edgeThreshold+1, nPairs*collisionProb(umiLen, edgeThreshold+1))
			}
		}
	}

	// Compute all-pairs edit distances; collect edges within the
	// (possibly adaptive) threshold. The Ukkonen cutoff bails out
	// after a few DP rows on dissimilar pairs, so the bound directly
	// affects performance — a tighter edgeThreshold means faster
	// edge-finding.
	var edges []umiEdge

	if numThreads <= 1 || n < 4 {
		var buf levBuf
		for i := range n {
			for j := i + 1; j < n; j++ {
				dist := distFn(normalized[i], normalized[j], &buf, edgeThreshold)
				if dist <= edgeThreshold {
					edges = append(edges, umiEdge{i, j, dist})
				}
			}
		}
	} else {
		// Each worker gets a slice of rows and collects edges locally.
		// Pre-allocated DP buffers per worker avoid heap churn.
		workerEdges := make([][]umiEdge, numThreads)
		var wg sync.WaitGroup

		// Distribute rows round-robin so work is roughly balanced
		// (earlier rows have more pairs than later ones).
		for w := 0; w < numThreads; w++ {
			wg.Add(1)
			go func(workerID int) {
				defer wg.Done()
				var buf levBuf
				var local []umiEdge
				for i := workerID; i < n; i += numThreads {
					for j := i + 1; j < n; j++ {
						dist := distFn(normalized[i], normalized[j], &buf, edgeThreshold)
						if dist <= edgeThreshold {
							local = append(local, umiEdge{i, j, dist})
						}
					}
				}
				workerEdges[workerID] = local
			}(w)
		}
		wg.Wait()

		for _, localEdges := range workerEdges {
			edges = append(edges, localEdges...)
		}
	}

	// ---------------------------------------------------------------
	// Cluster UMIs from edges using the selected method.
	//
	// All methods produce the same output: compMembers, a map from
	// cluster ID → list of member indices into umis/normalized.
	//
	//   connected  : single-linkage via union-find on all edges.
	//                Can chain (A↔B↔C even if d(A,C) >> threshold).
	//   adjacency  : greedy assignment. Highest-count UMI becomes a
	//                center, all its direct neighbors join. No chaining.
	//   directional: filter edges by PCR error count model, then
	//                union-find on filtered edges. Only low-count UMIs
	//                merge into high-count ones.
	//   tiered     : BFS from centers with decreasing edit distance
	//                per hop. More permissive than adjacency (multiple
	//                hops) but prevents chaining at high distances.
	// ---------------------------------------------------------------
	compMembers := make(map[int][]int)

	switch umiClusterMethod {
	case "adjacency":
		// Build adjacency list from edges.
		neighbors := make(map[int][]int, n)
		for _, e := range edges {
			neighbors[e.i] = append(neighbors[e.i], e.j)
			neighbors[e.j] = append(neighbors[e.j], e.i)
		}
		// Greedy assignment: process UMIs in count-descending order
		// (umis is already sorted by count desc, ties by umi asc).
		// Each unassigned UMI becomes a cluster center; all its
		// unassigned direct neighbors join that cluster.
		clusterOf := make([]int, n)
		for i := range clusterOf {
			clusterOf[i] = -1
		}
		for i := range umis {
			if clusterOf[i] != -1 {
				continue
			}
			clusterOf[i] = i // new cluster center
			for _, j := range neighbors[i] {
				if clusterOf[j] == -1 {
					clusterOf[j] = i
				}
			}
		}
		for i, center := range clusterOf {
			if center == -1 {
				center = i // isolated UMI, no edges
			}
			compMembers[center] = append(compMembers[center], i)
		}

	case "directional":
		// Filter edges by PCR error count model: only keep edge (i,j)
		// if the lower-count UMI could plausibly be a PCR/sequencing
		// error of the higher-count one. The formula from UMI-tools
		// (Smith et al. 2017) is:
		//   count(low) ≤ 2 × count(high) × (1/4)^distance
		// which models the expected number of errors at `distance`
		// substitutions with a per-base error rate of 1/4.
		var filtered []umiEdge
		for _, e := range edges {
			lo, hi := e.i, e.j
			if umis[lo].count > umis[hi].count {
				lo, hi = hi, lo
			}
			maxErrorCount := 2.0 * float64(umis[hi].count) * math.Pow(0.25, float64(e.dist))
			if float64(umis[lo].count) <= maxErrorCount {
				filtered = append(filtered, e)
			}
		}
		// Union-find on the filtered edges.
		parent := make([]int, n)
		for i := range parent {
			parent[i] = i
		}
		var find func(int) int
		find = func(x int) int {
			if parent[x] != x {
				parent[x] = find(parent[x])
			}
			return parent[x]
		}
		for _, e := range filtered {
			px, py := find(e.i), find(e.j)
			if px != py {
				parent[px] = py
			}
		}
		for i := range umis {
			compMembers[find(i)] = append(compMembers[find(i)], i)
		}

	case "tiered":
		// BFS from high-count centers with decreasing edit distance
		// threshold at each hop. With --umi-edit-distance T:
		//   Hop 0: center (highest-count unassigned UMI)
		//   Hop 1: neighbors at d ≤ T
		//   Hop 2: neighbors at d ≤ T-1
		//   Hop k: neighbors at d ≤ T-k (stop when T-k < 1)
		//
		// This is more permissive than adjacency (which is one hop)
		// but prevents chaining at high distances: you can't go
		// center→A(d=3)→B(d=3) because hop 2 requires d≤2.

		// Build adjacency list with distance annotations.
		type neighbor struct {
			idx  int
			dist int
		}
		neighbors := make(map[int][]neighbor, n)
		for _, e := range edges {
			neighbors[e.i] = append(neighbors[e.i], neighbor{e.j, e.dist})
			neighbors[e.j] = append(neighbors[e.j], neighbor{e.i, e.dist})
		}

		clusterOf := make([]int, n)
		for i := range clusterOf {
			clusterOf[i] = -1
		}

		// Process UMIs in count-descending order (umis already sorted).
		for i := range umis {
			if clusterOf[i] != -1 {
				continue
			}
			// This UMI is a new cluster center.
			clusterOf[i] = i

			// BFS with decreasing threshold per hop.
			type bfsEntry struct {
				idx int
				hop int
			}
			queue := []bfsEntry{{i, 0}}
			for len(queue) > 0 {
				cur := queue[0]
				queue = queue[1:]
				nextHop := cur.hop + 1
				maxDist := edgeThreshold - cur.hop
				if maxDist < 1 {
					continue // no more hops allowed
				}
				for _, nb := range neighbors[cur.idx] {
					if clusterOf[nb.idx] != -1 {
						continue
					}
					if nb.dist <= maxDist {
						clusterOf[nb.idx] = i
						queue = append(queue, bfsEntry{nb.idx, nextHop})
					}
				}
			}
		}

		for i, center := range clusterOf {
			if center == -1 {
				center = i
			}
			compMembers[center] = append(compMembers[center], i)
		}

	default: // "connected" — current single-linkage behavior
		parent := make([]int, n)
		for i := range parent {
			parent[i] = i
		}
		var find func(int) int
		find = func(x int) int {
			if parent[x] != x {
				parent[x] = find(parent[x])
			}
			return parent[x]
		}
		for _, e := range edges {
			px, py := find(e.i), find(e.j)
			if px != py {
				parent[px] = py
			}
		}
		for i := range umis {
			compMembers[find(i)] = append(compMembers[find(i)], i)
		}
	}

	// Compute max pairwise edit distance within each component.
	//
	// Single-linkage clustering means a component can be wider than the
	// per-edge threshold — chained A–B–C with d(A,B) <= T and d(B,C) <= T
	// can still have d(A,C) > T. We report that width, capped at
	// 3 x umiClusterEditThreshold to avoid unbounded DPs.
	//
	// For large clusters this is O(cluster²) — the same cost structure as
	// the edge-finding phase — so we parallelize it the same way:
	// round-robin rows across workers with per-worker DP buffers.
	compMaxDist := make(map[int]int)
	intraMaxBound := 3 * umiClusterEditThreshold
	for root, indices := range compMembers {
		nIdx := len(indices)
		if nIdx <= 1 {
			continue
		}
		if nIdx > 10000 {
			// Too expensive even parallelized — O(n²) on 10k+ members
			// would take tens of minutes. Report -1 so the caller knows
			// the value was skipped rather than actually zero.
			compMaxDist[root] = -1
			continue
		}
		if numThreads <= 1 || nIdx < 4 {
			// Small cluster or single-threaded: serial path.
			var buf levBuf
			maxDist := 0
			for a := 0; a < nIdx; a++ {
				for b := a + 1; b < nIdx; b++ {
					d := distFn(normalized[indices[a]], normalized[indices[b]], &buf, intraMaxBound)
					if d > maxDist {
						maxDist = d
					}
				}
			}
			compMaxDist[root] = maxDist
		} else {
			// Large cluster: distribute rows round-robin across workers.
			workerMax := make([]int, numThreads)
			var wg sync.WaitGroup
			for w := 0; w < numThreads; w++ {
				wg.Add(1)
				go func(workerID int) {
					defer wg.Done()
					var buf levBuf
					localMax := 0
					for a := workerID; a < nIdx; a += numThreads {
						for b := a + 1; b < nIdx; b++ {
							d := distFn(normalized[indices[a]], normalized[indices[b]], &buf, intraMaxBound)
							if d > localMax {
								localMax = d
							}
						}
					}
					workerMax[workerID] = localMax
				}(w)
			}
			wg.Wait()
			maxDist := 0
			for _, m := range workerMax {
				if m > maxDist {
					maxDist = m
				}
			}
			compMaxDist[root] = maxDist
		}
	}

	// Pick the representative for each cluster: the highest-count member.
	// Since umis is sorted by count desc (ties broken by umi asc), the
	// member with the smallest index is the highest-count UMI in the
	// cluster. We normalize its separator for consistent output.
	compRepresentative := make(map[int]string)
	for root, indices := range compMembers {
		bestIdx := indices[0]
		for _, idx := range indices[1:] {
			if idx < bestIdx {
				bestIdx = idx
			}
		}
		compRepresentative[root] = normalizeUMISeparator(umis[bestIdx].umi)
	}

	// Build results and populate globalRepresentative map. We iterate
	// compMembers (not a find() call) because the union-find is only
	// defined inside the switch branches above.
	results := make([]umiClusterResult, n)
	for root, indices := range compMembers {
		cons := compRepresentative[root]
		maxDist := compMaxDist[root]
		for _, k := range indices {
			globalRepresentative[umis[k].umi] = cons
			results[k] = umiClusterResult{
				umi:               umis[k].umi,
				representative:    cons,
				count:             umis[k].count,
				maxIntraClustDist: maxDist,
			}
		}
	}
	return results
}

// clusterStats returns the number of distinct representative UMIs
// (cluster count) and the size of the largest cluster (by number of
// unique UMI members) in a single pass over the results.
func clusterStats(results []umiClusterResult) (numClusters int, maxClusterSize int) {
	counts := make(map[string]int)
	for _, r := range results {
		counts[r.representative]++
	}
	for _, c := range counts {
		if c > maxClusterSize {
			maxClusterSize = c
		}
	}
	return len(counts), maxClusterSize
}

var umiClusterOutput string
var umiClusterTag string
var umiClusterOrigTag string
var umiClusterOverlap int
var umiClusterWholeGenome bool
var umiClusterNoStrand bool
var umiClusterEditThreshold int
var umiClusterCountsFilename string
var umiClusterMI bool
var umiClusterMatchOneEnd bool
var umiClusterThreads int
var umiClusterSkipRefs string
var umiClusterRegion string
var umiClusterHPDist bool
var umiClusterMethod string
var umiClusterAdaptiveThreshold bool
var umiClusterMaxFalsePairs float64

func init() {
	ontUmiClusterCmd.Flags().StringVarP(&umiClusterOutput, "output", "o", "", "Output BAM file path (required)")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterTag, "umi-tag", "RX", "SAM tag containing UMI sequence")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterOrigTag, "orig-umi-tag", "OX", "SAM tag to store original UMI before correction")
	ontUmiClusterCmd.Flags().IntVar(&umiClusterOverlap, "overlap", 50, "Maximum gap (bp) between reads to group them together")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterWholeGenome, "whole-genome", false, "Process all UMIs as a single group (ignore coordinates)")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterNoStrand, "no-strand", false, "Ignore strand when grouping reads (default: group by strand)")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterSkipRefs, "ignore-refs", "", "References to ignore (reads will be passed through with original UMI) (comma-separated)")
	ontUmiClusterCmd.Flags().IntVar(&umiClusterEditThreshold, "umi-edit-distance", 3, "Maximum Levenshtein edit distance to cluster two UMIs")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterCountsFilename, "umi-counts", "", "Write per-component UMI summary to this file")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterMI, "mi", false, "Add MI tag with molecule group ID to output reads")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterMatchOneEnd, "match-one-end", false, "Match reads if EITHER 5' or 3' ends are within gap (default: require BOTH ends)")
	ontUmiClusterCmd.Flags().IntVarP(&umiClusterThreads, "threads", "t", 1, "Threads for UMI clustering")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterRegion, "region", "", "Process only this region (e.g. 'chr19' or 'chr19:1000-2000'); disables the skipped-ref and unmapped passes")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterHPDist, "hp-dist", false, "Use HP-aware edit distance where homopolymer indels cost 0 (handles ONT HP-length errors without distorting non-HP mutations)")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterMethod, "umi-cluster-method", "adjacency", "UMI clustering method: connected (single-linkage), adjacency (greedy, no chaining), directional (PCR error count model), tiered (distance-attenuated BFS clustering)")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterAdaptiveThreshold, "adaptive-threshold", false, "Reduce edit distance threshold for large components to limit random collisions")
	ontUmiClusterCmd.Flags().Float64Var(&umiClusterMaxFalsePairs, "max-false-pairs", 1.0, "Maximum expected random false pairs allowed before reducing threshold (used with --adaptive-threshold)")
}
