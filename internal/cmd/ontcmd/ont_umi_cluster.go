package ontcmd

import (
	"fmt"
	"io"
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

		inputFile := args[0]

		skipRefs := []string{}
		if umiClusterSkipRefs != "" {
			skipRefs = strings.Split(umiClusterSkipRefs, ",")
		}

		// Open counts writer if requested
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

// groupWorkItem is sent to worker goroutines for parallel UMI clustering.
// It represents a read-overlap-group: reads with similar 5'/3' ends on the
// same strand. Within each group, UMIs are clustered by edit distance.
type groupWorkItem struct {
	recs     []*htsio.SamRecord
	rname    string
	strand   string
	minStart int
	maxEnd   int
}

// groupResult is returned by workers after UMI clustering.
type groupResult struct {
	item      groupWorkItem
	representative map[string]string
	results   []umiClusterResult
}

// umiClusterOverlapMode groups reads by 5' and/or 3' end proximity using a
// buffer + union-find approach, then clusters UMIs within each component.
// Chromosomes are processed in parallel (largest first) to maximize throughput.
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

	// Build skip set for fast lookup.
	skipSet := make(map[string]bool)
	for _, r := range skipRefs {
		skipSet[r] = true
	}

	// Sort chromosomes by length (largest first).
	refs := header.References()
	type chromJob struct {
		name   string
		length int
		idx    int // original header order
	}
	var jobs []chromJob
	for i, ref := range refs {
		if skipSet[ref.Name] {
			continue
		}
		jobs = append(jobs, chromJob{name: ref.Name, length: ref.Length, idx: i})
	}
	sort.Slice(jobs, func(i, j int) bool {
		return jobs[i].length > jobs[j].length
	})

	// Open writer — samtools sort handles merging unsorted chromosome output.
	writerThreads := 2
	if umiClusterThreads > 4 {
		writerThreads = 4
	}
	wopts := htsio.SamWriterOptions(header).BAM().SortCoord().Threads(writerThreads)
	writer, err := htsio.NewSamWriter(umiClusterOutput, wopts)
	if err != nil {
		return err
	}

	// Shared state.
	var nextMI atomic.Int64
	nextMI.Store(1)
	var totalReads atomic.Int64
	var totalChanged atomic.Int64
	var countsMu sync.Mutex // protects countsWriter

	if umiClusterThreads >= 4 {
		// Parallel: process threads-2 chromosomes concurrently.
		chromWorkers := umiClusterThreads - 2
		sem := make(chan struct{}, chromWorkers)

		var wg sync.WaitGroup
		var firstErr atomic.Value

		for _, job := range jobs {
			sem <- struct{}{} // acquire slot
			wg.Add(1)
			go func(chromName string) {
				defer wg.Done()
				defer func() { <-sem }() // release slot

				if err := processChromosome(inputFile, chromName, writer, &countsMu,
					&nextMI, &totalReads, &totalChanged, countsWriter); err != nil {
					firstErr.CompareAndSwap(nil, err)
				}
			}(job.name)
		}
		wg.Wait()

		if e := firstErr.Load(); e != nil {
			writer.Close()
			return e.(error)
		}
	} else {
		// Sequential: process one chromosome at a time.
		for _, job := range jobs {
			if err := processChromosome(inputFile, job.name, writer, &countsMu,
				&nextMI, &totalReads, &totalChanged, countsWriter); err != nil {
				writer.Close()
				return err
			}
		}
	}

	// Process skipped refs: pass through without UMI clustering.
	if len(skipRefs) > 0 {
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
				writer.Write(rec)
			}
			r.Close()
		}
	}

	// Process unmapped reads (unless --ignore-unmapped is set).
	if !umiClusterSkipUnmapped {
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
			writer.Write(rec)
		}
		unmappedReader.Close()
	}

	fmt.Fprintf(os.Stderr, "Total reads: %d, UMIs corrected: %d\n",
		totalReads.Load(), totalChanged.Load())
	return writer.Close()
}

// processChromosome runs the overlap-mode grouping + UMI clustering for a
// single chromosome. The writer is thread-safe; countsMu protects countsWriter.
func processChromosome(
	inputFile string,
	chromName string,
	writer *htsio.SamtoolsSamWriter,
	countsMu *sync.Mutex,
	nextMI *atomic.Int64,
	totalReads *atomic.Int64,
	totalChanged *atomic.Int64,
	countsWriter io.Writer,
) error {
	ropts := htsio.NewSamReaderOpts().Region(chromName)
	reader, err := htsio.NewSamReader(inputFile, ropts)
	if err != nil {
		return err
	}
	defer reader.Close()

	fmt.Fprintf(os.Stderr, "Processing %s...\n", chromName)

	// Per-chromosome state.
	var buffer []*bufferedRead
	uf := newUnionFind(1024)
	activeCount := make(map[int]int)
	componentReads := make(map[int][]*bufferedRead)
	globalID := 0

	// Collect results locally, then write in batch to reduce lock contention.
	type processedGroup struct {
		recs           []*htsio.SamRecord
		countsLines    []string
		readsChanged   int
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
		results := clusterUMIsParallel(umiCounts, representative, umiVerbose, 1)

		representativeCount := countRepresentative(results)
		strandLabel := ""
		if !umiClusterNoStrand {
			strandLabel = "(" + strand + ") "
		}
		fmt.Fprintf(os.Stderr, "%s:%d-%d: %s%d reads, %d unique UMIs -> %d representative\n",
			rname, minStart, maxEnd, strandLabel,
			len(recs), len(results), representativeCount)

		// Assign MI values using the shared atomic counter.
		repToMI := make(map[string]string)
		if umiClusterMI || countsWriter != nil {
			for _, r := range results {
				if _, ok := repToMI[r.representative]; !ok {
					mi := nextMI.Add(1) - 1
					repToMI[r.representative] = fmt.Sprintf("mi_%09d", mi)
				}
			}
		}

		// Build counts lines.
		item := &groupWorkItem{
			recs:     recs,
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
			recs:        recs,
			countsLines: cl,
			readsChanged: changed,
		}
	}

	writeGroup := func(pg *processedGroup) {
		if pg == nil {
			return
		}
		for _, rec := range pg.recs {
			writer.Write(rec)
		}
		if countsWriter != nil && len(pg.countsLines) > 0 {
			countsMu.Lock()
			for _, line := range pg.countsLines {
				fmt.Fprintln(countsWriter, line)
			}
			countsMu.Unlock()
		}
		totalReads.Add(int64(len(pg.recs)))
		totalChanged.Add(int64(pg.readsChanged))
	}

	submitComponent := func(root int) {
		reads := componentReads[root]
		delete(componentReads, root)
		delete(activeCount, root)
		pg := processGroup(reads)
		writeGroup(pg)
	}

	ejectRead := func(b *bufferedRead) {
		root := uf.find(b.id)
		componentReads[root] = append(componentReads[root], b)
		activeCount[root]--
		if activeCount[root] == 0 {
			submitComponent(root)
		}
	}

	ejectAll := func() {
		for _, b := range buffer {
			ejectRead(b)
		}
		buffer = buffer[:0]
	}

	ejectExpired := func(curStart int) {
		kept := 0
		for _, b := range buffer {
			shouldEject := false
			if umiClusterMatchOneEnd {
				shouldEject = curStart > b.end+umiClusterOverlap
			} else {
				shouldEject = curStart-b.start > umiClusterOverlap
			}
			if shouldEject {
				ejectRead(b)
			} else {
				buffer[kept] = b
				kept++
			}
		}
		for i := kept; i < len(buffer); i++ {
			buffer[i] = nil
		}
		buffer = buffer[:kept]
	}

	mergeComponents := func(newRoot, oldRoot int) {
		activeCount[newRoot] += activeCount[oldRoot]
		delete(activeCount, oldRoot)
		if reads, ok := componentReads[oldRoot]; ok {
			componentReads[newRoot] = append(componentReads[newRoot], reads...)
			delete(componentReads, oldRoot)
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
		rec, err := reader.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		if rec.IsUnmapped() || rec.Cigar == "*" {
			// Unmapped reads in a region query — pass through.
			writeGroup(&processedGroup{recs: []*htsio.SamRecord{rec}})
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

		for _, b := range buffer {
			if b.strand != br.strand {
				continue
			}

			var matches bool
			if umiClusterMatchOneEnd {
				fivePrime := br.start-b.start <= umiClusterOverlap
				diff := br.end - b.end
				if diff < 0 {
					diff = -diff
				}
				threePrime := diff <= umiClusterOverlap
				matches = fivePrime || threePrime
			} else {
				diff := br.end - b.end
				if diff < 0 {
					diff = -diff
				}
				matches = diff <= umiClusterOverlap
			}

			if matches {
				newRoot, oldRoot, merged := uf.union(br.id, b.id)
				if merged {
					mergeComponents(newRoot, oldRoot)
				}
			}
		}

		buffer = append(buffer, br)
	}

	ejectAll()

	fmt.Fprintf(os.Stderr, "Finished %s\n", chromName)
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

	for {
		rec, err := reader.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if umiClusterSkipUnmapped && (rec.IsUnmapped() || rec.Cigar == "*") {
			continue
		}
		for i := range skipRefs {
			if rec.RefName == skipRefs[i] {
				continue
			}
		}
		if umi := getUMI(rec); umi != "" {
			umiCounts[umi]++
		}
	}
	reader.Close()

	// Cluster
	globalRepresentative := make(map[string]string)
	results := clusterUMIs(umiCounts, globalRepresentative, umiVerbose)
	representativeCount := countRepresentative(results)
	fmt.Fprintf(os.Stderr, "whole-genome: %d unique UMIs -> %d representative\n", len(umiCounts), representativeCount)

	// Pass 2: rewrite BAM
	addUMIClusterPGLine(header)

	reader2, err := htsio.NewSamReader(inputFile)
	if err != nil {
		return err
	}

	opts := htsio.SamWriterOptions(header).BAM().SortCoord()

	if umiClusterThreads > 1 {
		opts = opts.Threads(2)
	}

	writer, err := htsio.NewSamWriter(umiClusterOutput, opts)
	if err != nil {
		return err
	}

	changed := 0
	total := 0
	for {
		rec, err := reader2.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		total++

		if umiClusterSkipUnmapped && (rec.IsUnmapped() || rec.Cigar == "*") {
			writer.Write(rec)
			continue
		}
		for i := range skipRefs {
			if rec.RefName == skipRefs[i] {
				writer.Write(rec)
				continue
			}
		}

		if updateRecordUMI(rec, globalRepresentative) {
			changed++
		}

		if err := writer.Write(rec); err != nil {
			return err
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

	var lines []string
	for _, ci := range clusters {
		lines = append(lines, fmt.Sprintf("%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%s",
			item.rname, item.minStart, item.maxEnd, ci.mi, item.strand,
			ci.representative, ci.numReads, len(ci.umis), ci.maxEditDist,
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

// levDist computes Levenshtein edit distance reusing the provided buffers.
func levDist(a, b string, buf *levBuf) int {
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
		for j := 1; j <= n; j++ {
			if a[i-1] == b[j-1] {
				buf.curr[j] = buf.prev[j-1]
			} else {
				buf.curr[j] = 1 + min(buf.prev[j], buf.curr[j-1], buf.prev[j-1])
			}
		}
		buf.prev, buf.curr = buf.curr, buf.prev
	}
	return buf.prev[n]
}

// computeRepresentativeUMI picks the representative UMI from a set of cluster
// members using multiple sequence alignment followed by majority vote per
// computeRepresentativeUMI picks the representative UMI for a cluster.
// The most common UMI (by read count) is chosen. Ties are broken by longer
// normalized length, then by lexicographic order.
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
	editDist          int // Levenshtein edit distance from this UMI to the cluster representative
	maxIntraClustDist int // maximum pairwise edit distance between any two members of this cluster
}

type umiCount struct {
	umi   string
	count int
}

// clusterUMIs clusters UMIs using the configured thread count for the
// all-pairs edit distance computation (suitable for whole-genome mode
// where there is a single large group).
func clusterUMIs(umiCounts map[string]int, globalRepresentative map[string]string, verbose bool) []umiClusterResult {
	return clusterUMIsParallel(umiCounts, globalRepresentative, verbose, umiClusterThreads)
}

// umiEdge represents a pair of UMIs within the edit distance threshold.
type umiEdge struct{ i, j, dist int }

func clusterUMIsParallel(umiCounts map[string]int, globalRepresentative map[string]string, verbose bool, numThreads int) []umiClusterResult {
	if len(umiCounts) <= 1 {
		// Single UMI or empty, nothing to cluster
		var results []umiClusterResult
		for umi, count := range umiCounts {
			globalRepresentative[umi] = umi
			results = append(results, umiClusterResult{umi: umi, representative: umi, count: count})
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
	var edges []umiEdge

	if numThreads <= 1 || n < 4 {
		var buf levBuf
		for i := range n {
			for j := i + 1; j < n; j++ {
				dist := levDist(normalized[i], normalized[j], &buf)
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
						dist := levDist(normalized[i], normalized[j], &buf)
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

	// Group UMIs by component root; track max intra-cluster distance per component.
	compMembers := make(map[int][]umiCount)
	for i := range umis {
		root := find(i)
		compMembers[root] = append(compMembers[root], umis[i])
	}
	// Compute max pairwise edit distance within each component.
	compMaxDist := make(map[int]int)
	for root, compMems := range compMembers {
		if len(compMems) <= 1 {
			continue
		}
		// Collect indices into the normalized array for this component.
		var indices []int
		for i := range umis {
			if find(i) == root {
				indices = append(indices, i)
			}
		}
		var buf levBuf
		maxDist := 0
		for a := 0; a < len(indices); a++ {
			for b := a + 1; b < len(indices); b++ {
				d := levDist(normalized[indices[a]], normalized[indices[b]], &buf)
				if d > maxDist {
					maxDist = d
				}
			}
		}
		compMaxDist[root] = maxDist
	}

	// Compute representative per component (most common UMI).
	compRepresentative := make(map[int]string)
	for root, members := range compMembers {
		compRepresentative[root] = computeRepresentativeUMI(members)
	}

	// Build results and populate globalRepresentative map.
	var buf levBuf
	results := make([]umiClusterResult, n)
	for k := range umis {
		root := find(k)
		cons := compRepresentative[root]
		globalRepresentative[umis[k].umi] = cons
		results[k] = umiClusterResult{
			umi:               umis[k].umi,
			representative:    cons,
			count:             umis[k].count,
			editDist:          levDist(normalized[k], normalizeUMISeparator(cons), &buf),
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
var umiClusterSkipUnmapped bool

var umiVerbose bool

func init() {
	ontUmiClusterCmd.Flags().StringVarP(&umiClusterOutput, "output", "o", "", "Output BAM file path (required)")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterTag, "umi-tag", "RX", "SAM tag containing UMI sequence")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterOrigTag, "orig-umi-tag", "OX", "SAM tag to store original UMI before correction")
	ontUmiClusterCmd.Flags().IntVar(&umiClusterOverlap, "overlap", 50, "Maximum gap (bp) between reads to group them together")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterWholeGenome, "whole-genome", false, "Process all UMIs as a single group (ignore coordinates)")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterNoStrand, "no-strand", false, "Ignore strand when grouping reads (default: group by strand)")
	ontUmiClusterCmd.Flags().BoolVarP(&umiVerbose, "verbose", "v", false, "Verbose logging")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterSkipRefs, "ignore-refs", "", "References to ignore (reads will be passed through with original UMI) (comma-separated)")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterSkipUnmapped, "ignore-unmapped", false, "Ignore unmapped reads (reads will be passed through with original UMI)")
	ontUmiClusterCmd.Flags().IntVar(&umiClusterEditThreshold, "umi-edit-distance", 3, "Maximum Levenshtein edit distance to cluster two UMIs")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterCountsFilename, "umi-counts", "", "Write per-component UMI summary to this file")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterMI, "mi", false, "Add MI tag with molecule group ID to output reads")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterMatchOneEnd, "match-one-end", false, "Match reads if EITHER 5' or 3' ends are within gap (default: require BOTH ends)")
	ontUmiClusterCmd.Flags().IntVarP(&umiClusterThreads, "threads", "t", 1, "Threads for UMI clustering")
}
