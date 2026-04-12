package ontcmd

import (
	"fmt"
	"io"
	"os"
	"sort"
	"strings"
	"sync"

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
			if _, err := fmt.Fprintln(countsWriter, umiCountsHeaderLine); err != nil {
				return fmt.Errorf("writing umi-counts header: %w", err)
			}
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

// groupWorkItem describes a single read-overlap-group: reads sharing
// 5'/3' endpoints within the overlap threshold on the same strand. Used
// to carry the coordinate-level metadata needed by buildUMIClusterCounts
// (chrom, strand, span); the underlying SamRecords are passed separately
// and aren't stored here.
type groupWorkItem struct {
	rname    string
	strand   string
	minStart int
	maxEnd   int
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
	reader, err := htsio.NewSamReader(inputFile)
	if err != nil {
		return err
	}
	header, err := reader.Header()
	if err != nil {
		return fmt.Errorf("failed to read header: %w", err)
	}
	reader.Close()
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

	// Per-run state. Plain ints/counters are safe because everything below
	// this point runs on the caller's single goroutine — we removed the
	// per-chromosome parallelism in favor of one-region-per-invocation.
	// If a future rework reintroduces parallel regions, these will need
	// to become atomics (or the counters need to move into a synchronized
	// aggregator), but for now the simpler types avoid pointless atomic
	// ops on every record written.
	var nextMI int64 = 1
	var totalReads int64
	var totalChanged int64

	// Process each region sequentially. processRegion handles the buffer
	// + union-find grouping and per-group UMI clustering with the full
	// --threads budget for clusterUMIs.
	for _, region := range regions {
		if err := processRegion(inputFile, region, writer,
			&nextMI, &totalReads, &totalChanged, countsWriter); err != nil {
			writer.Close()
			return err
		}
	}

	// In full-file mode (no --region), pass through skipped refs and
	// unmapped reads unchanged so the output BAM is complete. In
	// --region mode, the user is orchestrating per-region jobs
	// externally, so these pass-through passes are someone else's
	// concern — we skip them here to avoid duplicating records.
	if umiClusterRegion == "" {
		// Skipped refs: pass through without UMI clustering.
		for _, refName := range skipRefs {
			ropts := htsio.NewSamReaderOpts().Region(refName)
			r, err := htsio.NewSamReader(inputFile, ropts)
			if err != nil {
				writer.Close()
				return err
			}
			for {
				rec, err := r.Next()
				if err == io.EOF {
					break
				}
				if err != nil {
					r.Close()
					writer.Close()
					return err
				}
				if err := writer.Write(rec); err != nil {
					r.Close()
					writer.Close()
					return fmt.Errorf("writing skipped-ref record: %w", err)
				}
			}
			r.Close()
		}

		// Unmapped reads: pass through unchanged. They have no
		// position so there's nothing to cluster against, but they
		// belong in the output BAM so downstream tools see a
		// complete file.
		unmappedOpts := htsio.NewSamReaderOpts().FlagRequired(0x4)
		unmappedReader, err := htsio.NewSamReader(inputFile, unmappedOpts)
		if err != nil {
			writer.Close()
			return err
		}
		for {
			rec, err := unmappedReader.Next()
			if err == io.EOF {
				break
			}
			if err != nil {
				unmappedReader.Close()
				writer.Close()
				return err
			}
			if err := writer.Write(rec); err != nil {
				unmappedReader.Close()
				writer.Close()
				return fmt.Errorf("writing unmapped record: %w", err)
			}
		}
		unmappedReader.Close()
	}

	fmt.Fprintf(os.Stderr, "Total reads: %d, UMIs corrected: %d\n",
		totalReads, totalChanged)
	return writer.Close()
}

// processRegion runs the overlap-mode grouping + UMI clustering for a
// single samtools region (a chromosome name like "chr19" or a region like
// "chr19:1000000-2000000"). The caller (umiClusterOverlapMode) invokes
// this sequentially per region, so the counter pointers are plain *int64
// rather than atomic — there's no concurrent mutation.
func processRegion(
	inputFile string,
	region string,
	writer *htsio.SamtoolsSamWriter,
	nextMI *int64,
	totalReads *int64,
	totalChanged *int64,
	countsWriter io.Writer,
) error {
	ropts := htsio.NewSamReaderOpts().Region(region)
	reader, err := htsio.NewSamReader(inputFile, ropts)
	if err != nil {
		return err
	}
	defer reader.Close()

	fmt.Fprintf(os.Stderr, "Processing %s...\n", region)

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

	// processedGroup holds a fully-clustered read-overlap-group: the
	// final SamRecords to write, the counts-file lines that describe
	// this group's clusters, and the number of records whose UMI tag
	// was rewritten. The split between processGroup (which produces
	// one of these) and writeGroup (which consumes it) exists so the
	// per-group stderr log line appears before the records hit the
	// writer, which makes killed jobs easier to diagnose.
	type processedGroup struct {
		recs         []*htsio.SamRecord
		countsLines  []string
		readsChanged int
	}

	processGroup := func(reads []*bufferedRead) *processedGroup {
		if len(reads) == 0 {
			return nil
		}
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

		// Cluster UMIs.
		umiCounts := make(map[string]int)
		for _, rec := range recs {
			if umi := getUMI(rec); umi != "" {
				umiCounts[umi]++
			}
		}
		representative := make(map[string]string)
		// Use the full --threads budget for per-group clustering.
		results := clusterUMIs(umiCounts, representative, umiClusterThreads)

		representativeCount := countRepresentative(results)
		strandLabel := ""
		if !umiClusterNoStrand {
			strandLabel = "(" + strand + ") "
		}
		fmt.Fprintf(os.Stderr, "%s:%d-%d: %s%d reads, %d unique UMIs -> %d representative\n",
			rname, minStart, maxEnd, strandLabel,
			len(recs), len(results), representativeCount)

		// Assign MI values.
		repToMI := make(map[string]string)
		if umiClusterMI || countsWriter != nil {
			for _, r := range results {
				if _, ok := repToMI[r.representative]; !ok {
					mi := *nextMI
					*nextMI++
					repToMI[r.representative] = fmt.Sprintf("mi_%09d", mi)
				}
			}
		}

		// Build counts lines.
		item := &groupWorkItem{
			rname:    rname,
			strand:   strand,
			minStart: minStart,
			maxEnd:   maxEnd,
		}
		var cl []string
		if countsWriter != nil && len(results) > 0 {
			cl = buildUMIClusterCounts(item, results, repToMI)
		}

		// Update records.
		changed := 0
		for _, rec := range recs {
			origUMI := getUMI(rec)
			if updateRecordUMI(rec, representative) {
				changed++
			}
			if umiClusterMI && origUMI != "" {
				if rep, ok := representative[origUMI]; ok {
					if mi, ok := repToMI[rep]; ok {
						rec.Tags["MI"] = htsio.SamTag{Type: 'Z', Value: mi}
					}
				}
			}
		}

		return &processedGroup{
			recs:         recs,
			countsLines:  cl,
			readsChanged: changed,
		}
	}

	writeGroup := func(pg *processedGroup) error {
		if pg == nil {
			return nil
		}
		for _, rec := range pg.recs {
			if err := writer.Write(rec); err != nil {
				return fmt.Errorf("writing record to output BAM: %w", err)
			}
		}
		if countsWriter != nil && len(pg.countsLines) > 0 {
			for _, line := range pg.countsLines {
				if _, err := fmt.Fprintln(countsWriter, line); err != nil {
					return fmt.Errorf("writing umi-counts line: %w", err)
				}
			}
		}
		*totalReads += int64(len(pg.recs))
		*totalChanged += int64(pg.readsChanged)
		return nil
	}

	// writeErr is set by any writeGroup/submitComponent failure. The
	// closures check it (and no-op once set) and the main read loop
	// bails out on the next iteration.
	var writeErr error

	submitComponent := func(root int) {
		reads := componentReads[root]
		delete(componentReads, root)
		delete(activeCount, root)
		pg := processGroup(reads)
		if writeErr != nil {
			return
		}
		if err := writeGroup(pg); err != nil {
			writeErr = err
		}
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

	// ejectAll drains every remaining active read. Called after the last
	// record in the region has been read.
	ejectAll := func() {
		for _, b := range active {
			if writeErr != nil {
				return
			}
			// Inline the union-find bookkeeping (don't call
			// removeFromIndex — we're about to wipe the maps).
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

	for {
		if writeErr != nil {
			return writeErr
		}
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

		readStart := rec.Pos - 1
		readEnd := readStart + htsio.CigarRefLen(rec.Cigar)
		strand := readStrand(rec)

		ejectExpired(readStart)

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

		// Query the bin indices for overlapping reads. In default mode
		// (match-both-ends), all active reads already pass the 5' check
		// (the start-based ejection guarantees it), so we only query the
		// end index. In match-one-end mode, we query both indices and
		// take the union (via union-find short-circuit on duplicates).
		//
		// Each index query touches at most 3 bins (target ± 1), so the
		// work per read is O(matches) instead of O(buffer_size).

		endLo := (br.end - overlap) / overlap
		endHi := (br.end + overlap) / overlap

		if umiClusterMatchOneEnd {
			// 3'-end matches.
			for bin := endLo; bin <= endHi; bin++ {
				for _, b := range endIndex[bin] {
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
						}
					}
				}
			}
			// 5'-start matches: reads with start in
			// [br.start - overlap, br.start]. Since BAM is
			// coordinate-sorted, b.start <= br.start for all
			// buffered reads, so we query
			// [br.start - overlap, br.start].
			startLo := (br.start - overlap) / overlap
			startHi := br.start / overlap
			for bin := startLo; bin <= startHi; bin++ {
				for _, b := range startIndex[bin] {
					if b.strand != br.strand {
						continue
					}
					if br.start-b.start <= overlap {
						newRoot, oldRoot, merged := uf.union(br.id, b.id)
						if merged {
							mergeComponents(newRoot, oldRoot)
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
	if writeErr != nil {
		return writeErr
	}

	fmt.Fprintf(os.Stderr, "Finished %s\n", region)
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
	representativeCount := countRepresentative(results)
	fmt.Fprintf(os.Stderr, "whole-genome: %d unique UMIs -> %d representative\n", len(umiCounts), representativeCount)

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

// buildUMIClusterCounts returns tab-delimited counts lines for one read-overlap-group.
func buildUMIClusterCounts(item *groupWorkItem, results []umiClusterResult, repToMI map[string]string) []string {
	// Group results by representative UMI.
	type clusterInfo struct {
		mi             string
		representative string
		numReads       int      // total reads across all UMIs in this cluster
		umis           []string // all distinct original UMIs in this cluster
		maxEditDist    int
	}
	clusters := make(map[string]*clusterInfo)
	for _, r := range results {
		ci, ok := clusters[r.representative]
		if !ok {
			ci = &clusterInfo{
				mi:             repToMI[r.representative],
				representative: r.representative,
			}
			clusters[r.representative] = ci
		}
		ci.numReads += r.count
		ci.umis = append(ci.umis, r.umi)
		if r.maxIntraClustDist > ci.maxEditDist {
			ci.maxEditDist = r.maxIntraClustDist
		}
	}

	// Output format is BED6+ (standard BED6 columns followed by extras):
	//   chrom, start, end, name, score, strand,
	//   representative, numUMIs, maxEditDist, umis
	// where `name` is the molecule ID and `score` is the total read count
	// for the cluster. BED convention caps score at 1000, but downstream
	// tools (bedtools, IGV) accept larger values, and preserving the full
	// count is more useful than clamping.
	var lines []string
	for _, ci := range clusters {
		lines = append(lines, fmt.Sprintf("%s\t%d\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%s",
			item.rname, item.minStart, item.maxEnd, ci.mi, ci.numReads, item.strand,
			ci.representative, len(ci.umis), ci.maxEditDist,
			strings.Join(ci.umis, ",")))
	}
	return lines
}

// umiCountsHeaderLine is the commented column header written to the top of
// the --umi-counts file. Prefixed with '#' so BED parsers treat it as a
// comment and skip it.
const umiCountsHeaderLine = "#chrom\tstart\tend\tname\tscore\tstrand\trepresentative\tnumUMIs\tmaxEditDist\tumis"

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

	// Compute all-pairs edit distances; collect edges within threshold.
	// Partition rows across workers so each worker processes a contiguous
	// block of rows without channel overhead per pair.
	//
	// Every levDist call here is bounded by umiClusterEditThreshold: we
	// only care about pairs that pass the cluster threshold, so the
	// Ukkonen cutoff lets the DP bail out after just a few rows on
	// dissimilar UMIs. This is the main speedup for large components.
	var edges []umiEdge

	if numThreads <= 1 || n < 4 {
		var buf levBuf
		for i := range n {
			for j := i + 1; j < n; j++ {
				dist := levDist(normalized[i], normalized[j], &buf, umiClusterEditThreshold)
				if dist <= umiClusterEditThreshold {
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
						dist := levDist(normalized[i], normalized[j], &buf, umiClusterEditThreshold)
						if dist <= umiClusterEditThreshold {
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

	// Union-Find: find connected components.
	parent := make([]int, n)
	for i := range parent {
		parent[i] = i
	}
	var find func(int) int
	find = func(x int) int {
		if parent[x] != x {
			parent[x] = find(parent[x]) // path compression
		}
		return parent[x]
	}
	union := func(x, y int) {
		px, py := find(x), find(y)
		if px != py {
			parent[px] = py
		}
	}
	for _, e := range edges {
		union(e.i, e.j)
	}

	// Group UMIs by component root. We store indices into the `umis` /
	// `normalized` slices so we can look up both the umiCount and the
	// normalized string for a member without any additional passes.
	compMembers := make(map[int][]int)
	for i := range umis {
		root := find(i)
		compMembers[root] = append(compMembers[root], i)
	}

	// Compute max pairwise edit distance within each component.
	//
	// Single-linkage clustering means a component can be wider than the
	// per-edge threshold — chained A–B–C with d(A,B) <= T and d(B,C) <= T
	// can still have d(A,C) > T. We want to report that width, but we
	// cap it at 3 x umiClusterEditThreshold: beyond that we don't care
	// about the exact number, and capping keeps a pathological cluster
	// from wasting hours on unbounded DPs in a 100k-read component.
	//
	// The returned "max" is min(true_max, 3*threshold + 1) — callers can
	// treat any value of 3*threshold+1 as "at least 3*threshold+1".
	compMaxDist := make(map[int]int)
	intraMaxBound := 3 * umiClusterEditThreshold
	var buf levBuf
	for root, indices := range compMembers {
		if len(indices) <= 1 {
			continue
		}
		maxDist := 0
		for a := 0; a < len(indices); a++ {
			for b := a + 1; b < len(indices); b++ {
				d := levDist(normalized[indices[a]], normalized[indices[b]], &buf, intraMaxBound)
				if d > maxDist {
					maxDist = d
				}
			}
		}
		compMaxDist[root] = maxDist
	}

	// Compute representative per component (most common UMI). We build a
	// temporary umiCount slice from the stored indices; allocation is
	// negligible compared to the O(n^2) distance work above.
	compRepresentative := make(map[int]string)
	for root, indices := range compMembers {
		members := make([]umiCount, len(indices))
		for i, idx := range indices {
			members[i] = umis[idx]
		}
		compRepresentative[root] = computeRepresentativeUMI(members)
	}

	// Build results and populate globalRepresentative map.
	results := make([]umiClusterResult, n)
	for k := range umis {
		root := find(k)
		cons := compRepresentative[root]
		globalRepresentative[umis[k].umi] = cons
		results[k] = umiClusterResult{
			umi:               umis[k].umi,
			representative:    cons,
			count:             umis[k].count,
			maxIntraClustDist: compMaxDist[root],
		}
	}
	return results
}

func countRepresentative(results []umiClusterResult) int {
	seen := make(map[string]bool)
	for _, r := range results {
		seen[r.representative] = true
	}
	return len(seen)
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
}
