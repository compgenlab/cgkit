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

		inputFile := args[0]

		// Open counts writer if requested
		var countsWriter io.Writer
		var closeCounts func() error
		if umiClusterCountsFilename != "" {
			var err error
			countsWriter, closeCounts, err = openWriter(umiClusterCountsFilename)
			if err != nil {
				return fmt.Errorf("opening umi-counts: %w", err)
			}
			defer func() {
				if err := closeCounts(); err != nil {
					fmt.Fprintf(os.Stderr, "error closing umi-counts: %v\n", err)
				}
			}()
		}

		// Open BED writer if requested
		var bedWriter io.Writer
		var closeBed func() error
		if umiClusterBedFilename != "" {
			var err error
			bedWriter, closeBed, err = openWriter(umiClusterBedFilename)
			if err != nil {
				return fmt.Errorf("opening bed: %w", err)
			}
			defer func() {
				if err := closeBed(); err != nil {
					fmt.Fprintf(os.Stderr, "error closing bed: %v\n", err)
				}
			}()
		}

		if umiClusterWholeGenome {
			return umiClusterWholeGenomeMode(inputFile, countsWriter)
		}
		return umiClusterOverlapMode(inputFile, countsWriter, bedWriter)
	},
}

// overlapGroup tracks a set of buffered reads that overlap or are within
// the gap tolerance of each other, on the same strand (unless --no-strand).
type overlapGroup struct {
	rname    string
	strand   string // "+" or "-", or "" for --no-strand
	minStart int    // 0-based
	maxEnd   int    // 0-based, exclusive
	recs     []*htsio.SamRecord
}

func (g *overlapGroup) reset() {
	g.recs = nil
}

// readsNearby returns true if a read at readStart is within gap tolerance of this group.
// Reads that overlap OR whose gap to the group is <= tolerance are merged.
func (g *overlapGroup) readsNearby(rname string, strand string, readStart int) bool {
	if len(g.recs) == 0 {
		return false
	}
	if rname != g.rname {
		return false
	}
	if !umiClusterNoStrand && strand != g.strand {
		return false
	}
	// readStart is guaranteed >= g.minStart (coordinate sorted).
	// Gap = readStart - g.maxEnd (negative means overlap).
	// Merge if gap <= tolerance.
	return readStart <= g.maxEnd+umiClusterOverlap
}

// groupWorkItem is sent to worker goroutines for parallel UMI clustering.
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
	consensus map[string]string
	results   []umiClusterResult
}

// umiClusterOverlapMode processes each overlap group independently:
// buffer the group's reads, send to workers for UMI clustering in parallel,
// write results in order.
func umiClusterOverlapMode(inputFile string, countsWriter io.Writer, bedWriter io.Writer) error {
	reader, err := htsio.NewSamReader(inputFile)
	if err != nil {
		return err
	}
	if umiClusterThreads > 1 {
		reader.Threads(2)
	}

	// Read first record to populate header
	firstRec, err := reader.Next()
	if err != nil && err != io.EOF {
		return err
	}

	header := reader.Header()
	if header == nil {
		return fmt.Errorf("no header found in BAM file")
	}
	if err := validateCoordinateSorted(header); err != nil {
		return err
	}

	header.AddLine(fmt.Sprintf("@PG\tID:ont-umi-cluster\tPN:cgltk\tCL:ont-umi-cluster\tDS:UMI collapsing; consensus UMI written to %s, original preserved in %s", umiClusterTag, umiClusterOrigTag))

	writer, err := htsio.NewSamWriter(umiClusterOutput, header)
	if err != nil {
		return err
	}
	writer.Format(htsio.FormatBAM)
	if umiClusterThreads > 1 {
		writer.Threads(2)
	}

	numThreads := umiClusterThreads
	if numThreads <= 0 {
		numThreads = 1
	}

	type workItemWithCh struct {
		item     groupWorkItem
		resultCh chan *groupResult
	}

	workCh := make(chan workItemWithCh, numThreads)
	orderCh := make(chan chan *groupResult, numThreads*2)

	var workerWg sync.WaitGroup
	for range numThreads {
		workerWg.Add(1)
		go func() {
			defer workerWg.Done()
			for w := range workCh {
				umiCounts := make(map[string]int)
				for _, rec := range w.item.recs {
					if umi := getUMI(rec); umi != "" {
						umiCounts[umi]++
					}
				}
				consensus := make(map[string]string)
				results := clusterUMIs(umiCounts, consensus, umiVerbose)
				w.resultCh <- &groupResult{
					item:      w.item,
					consensus: consensus,
					results:   results,
				}
			}
		}()
	}

	// Writer goroutine: processes results in submission order.
	writerErrCh := make(chan error, 1)
	go func() {
		var writeErr error
		totalReads := 0
		totalChanged := 0
		regionCount := 0

		for resultCh := range orderCh {
			gr := <-resultCh
			if writeErr != nil {
				continue // drain but skip writing
			}

			if gr == nil {
				// This shouldn't happen in normal flow
				continue
			}

			// Write counts and BED
			regionCount++
			regionName := fmt.Sprintf("region_%d", regionCount)
			writeGroupCounts(countsWriter, gr.results, gr.item.rname, gr.item.minStart, gr.item.maxEnd, gr.item.strand, false)
			if bedWriter != nil {
				fmt.Fprintf(bedWriter, "%s\t%d\t%d\t%s\t0\t%s\n", gr.item.rname, gr.item.minStart, gr.item.maxEnd, regionName, gr.item.strand)
			}

			consensusCount := countConsensus(gr.results)
			strandLabel := ""
			if !umiClusterNoStrand {
				strandLabel = "(" + gr.item.strand + ") "
			}
			fmt.Fprintf(os.Stderr, "%s:%d-%d: %s%d reads, %d unique UMIs -> %d consensus\n",
				gr.item.rname, gr.item.minStart, gr.item.maxEnd, strandLabel,
				len(gr.item.recs), len(gr.results), consensusCount)

			// Update tags and write reads
			for _, rec := range gr.item.recs {
				tag, hasTag := rec.Tags[umiClusterTag]
				if hasTag {
					umi := tag.Value
					if cons, ok := gr.consensus[umi]; ok && cons != umi {
						rec.Tags[umiClusterOrigTag] = htsio.SamTag{Type: 'Z', Value: umi}
						rec.Tags[umiClusterTag] = htsio.SamTag{Type: 'Z', Value: cons}
						totalChanged++
					}
				}
				totalReads++
				if err := writer.Write(rec); err != nil {
					writeErr = err
				}
			}
		}

		fmt.Fprintf(os.Stderr, "Total reads: %d, UMIs corrected: %d\n", totalReads, totalChanged)
		writerErrCh <- writeErr
	}()

	// Maintain separate groups per strand.
	var plusGroup, minusGroup overlapGroup

	submitGroup := func(group *overlapGroup) {
		if len(group.recs) == 0 {
			return
		}
		item := groupWorkItem{
			recs:     group.recs,
			rname:    group.rname,
			strand:   group.strand,
			minStart: group.minStart,
			maxEnd:   group.maxEnd,
		}
		resultCh := make(chan *groupResult, 1)
		orderCh <- resultCh
		workCh <- workItemWithCh{item: item, resultCh: resultCh}
		group.recs = nil // detach slice so new group gets fresh storage
	}

	flushIfPast := func(group *overlapGroup, rname string, readStart int) {
		if len(group.recs) == 0 {
			return
		}
		if rname != group.rname || readStart > group.maxEnd+umiClusterOverlap {
			submitGroup(group)
		}
	}

	readStrand := func(rec *htsio.SamRecord) string {
		if rec.IsReverse() {
			return "-"
		}
		return "+"
	}

	addRecord := func(rec *htsio.SamRecord) {
		if rec.IsUnmapped() || rec.Cigar == "*" {
			// Unmapped reads: submit current groups first to maintain order,
			// then write through via a passthrough result.
			submitGroup(&plusGroup)
			if !umiClusterNoStrand {
				submitGroup(&minusGroup)
			}
			// Write unmapped read in order via a direct result
			resultCh := make(chan *groupResult, 1)
			orderCh <- resultCh
			resultCh <- &groupResult{
				item: groupWorkItem{recs: []*htsio.SamRecord{rec}},
			}
			return
		}

		readStart := rec.Pos - 1
		readEnd := readStart + htsio.CigarRefLen(rec.Cigar)
		strand := readStrand(rec)

		flushIfPast(&plusGroup, rec.RName, readStart)
		if !umiClusterNoStrand {
			flushIfPast(&minusGroup, rec.RName, readStart)
		}

		var group *overlapGroup
		if umiClusterNoStrand {
			group = &plusGroup
		} else if strand == "+" {
			group = &plusGroup
		} else {
			group = &minusGroup
		}

		if group.readsNearby(rec.RName, strand, readStart) {
			if readEnd > group.maxEnd {
				group.maxEnd = readEnd
			}
			group.recs = append(group.recs, rec)
		} else {
			submitGroup(group)
			group.rname = rec.RName
			if umiClusterNoStrand {
				group.strand = "."
			} else {
				group.strand = strand
			}
			group.minStart = readStart
			group.maxEnd = readEnd
			group.recs = append(group.recs, rec)
		}
	}

	if firstRec != nil {
		addRecord(firstRec)
	}
	for {
		rec, err := reader.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		addRecord(rec)
	}
	reader.Close()

	// Flush remaining groups
	submitGroup(&plusGroup)
	submitGroup(&minusGroup)

	close(workCh)
	workerWg.Wait()
	close(orderCh)

	writeErr := <-writerErrCh
	if writeErr != nil {
		return writeErr
	}
	return writer.Close()
}

// umiClusterWholeGenomeMode uses two passes over the entire file:
// pass 1 collects all UMI counts, pass 2 rewrites with consensus UMIs.
func umiClusterWholeGenomeMode(inputFile string, countsWriter io.Writer) error {
	// Pass 1: collect all UMIs
	reader, err := htsio.NewSamReader(inputFile)
	if err != nil {
		return err
	}
	umiCounts := make(map[string]int)

	firstRec, err := reader.Next()
	if err != nil && err != io.EOF {
		return err
	}
	header := reader.Header()
	if header == nil {
		return fmt.Errorf("no header found in BAM file")
	}

	if firstRec != nil {
		if umi := getUMI(firstRec); umi != "" {
			umiCounts[umi]++
		}
	}
	for {
		rec, err := reader.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if umi := getUMI(rec); umi != "" {
			umiCounts[umi]++
		}
	}
	reader.Close()

	// Cluster
	globalConsensus := make(map[string]string)
	results := clusterUMIs(umiCounts, globalConsensus, umiVerbose)
	writeGroupCounts(countsWriter, results, "", 0, 0, "", true)
	consensusCount := countConsensus(results)
	fmt.Fprintf(os.Stderr, "whole-genome: %d unique UMIs -> %d consensus\n", len(umiCounts), consensusCount)

	// Pass 2: rewrite BAM
	header.AddLine(fmt.Sprintf("@PG\tID:ont-umi-cluster\tPN:cgltk\tCL:ont-umi-cluster\tDS:UMI collapsing; consensus UMI written to %s, original preserved in %s", umiClusterTag, umiClusterOrigTag))

	reader2, err := htsio.NewSamReader(inputFile)
	if err != nil {
		return err
	}
	writer, err := htsio.NewSamWriter(umiClusterOutput, header)
	if err != nil {
		return err
	}
	writer.Format(htsio.FormatBAM)

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

		tag, hasTag := rec.Tags[umiClusterTag]
		if hasTag {
			umi := tag.Value
			if cons, ok := globalConsensus[umi]; ok && cons != umi {
				rec.Tags[umiClusterOrigTag] = htsio.SamTag{Type: 'Z', Value: umi}
				rec.Tags[umiClusterTag] = htsio.SamTag{Type: 'Z', Value: cons}
				changed++
			}
		}

		if err := writer.Write(rec); err != nil {
			return err
		}
	}
	reader2.Close()

	fmt.Fprintf(os.Stderr, "Total reads: %d, UMIs corrected: %d\n", total, changed)
	return writer.Close()
}

func writeGroupCounts(w io.Writer, results []umiClusterResult, rname string, start0, end0 int, strand string, wholeGenome bool) {
	if w == nil || len(results) == 0 {
		return
	}

	// Sum counts per consensus UMI
	consensusTotals := make(map[string]int)
	for _, r := range results {
		consensusTotals[r.consensus] += r.count
	}

	for _, r := range results {
		if !wholeGenome {
			fmt.Fprintf(w, "%s\t%d\t%d\t%s\t", rname, start0, end0, strand)
		}
		fmt.Fprintf(w, "%s\t%s\t%d\t%d\t%d\t%d\n", r.umi, r.consensus, r.count, consensusTotals[r.consensus], r.editDist, r.maxIntraClustDist)
	}
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

// umiLevenshtein computes the Levenshtein edit distance between two
// separator-normalized UMI strings.
func umiLevenshtein(a, b string) int {
	m, n := len(a), len(b)
	if m == 0 {
		return n
	}
	if n == 0 {
		return m
	}
	prev := make([]int, n+1)
	curr := make([]int, n+1)
	for j := 0; j <= n; j++ {
		prev[j] = j
	}
	for i := 1; i <= m; i++ {
		curr[0] = i
		for j := 1; j <= n; j++ {
			if a[i-1] == b[j-1] {
				curr[j] = prev[j-1]
			} else {
				curr[j] = 1 + min(prev[j], curr[j-1], prev[j-1])
			}
		}
		prev, curr = curr, prev
	}
	return prev[n]
}

// computeConsensusUMI calculates a consensus UMI string from a set of cluster
// members using multiple sequence alignment followed by majority vote per
// computeConsensusUMI picks the representative UMI for a cluster.
// The most common UMI (by read count) is chosen. Ties are broken by longer
// normalized length, then by lexicographic order.
func computeConsensusUMI(members []umiCount) string {
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
	consensus         string
	count             int
	editDist          int // Levenshtein edit distance from this UMI to the cluster consensus
	maxIntraClustDist int // maximum pairwise edit distance between any two members of this cluster
}

type umiCount struct {
	umi   string
	count int
}

func clusterUMIs(umiCounts map[string]int, globalConsensus map[string]string, verbose bool) []umiClusterResult {
	return clusterUMIsParallel(umiCounts, globalConsensus, verbose, umiClusterThreads)
}

// umiEdge represents a pair of UMIs within the edit distance threshold.
type umiEdge struct{ i, j, dist int }

func clusterUMIsParallel(umiCounts map[string]int, globalConsensus map[string]string, verbose bool, numThreads int) []umiClusterResult {
	if len(umiCounts) <= 1 {
		// Single UMI or empty, nothing to cluster
		var results []umiClusterResult
		for umi, count := range umiCounts {
			globalConsensus[umi] = umi
			results = append(results, umiClusterResult{umi: umi, consensus: umi, count: count})
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
	type pairJob struct{ i, j int }
	var edges []umiEdge

	if numThreads <= 1 {
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				dist := umiLevenshtein(normalized[i], normalized[j])
				if dist <= umiClusterEditThreshold {
					edges = append(edges, umiEdge{i, j, dist})
					if verbose {
						fmt.Printf("EDGE edit_dist=%d: %s -- %s\n", dist, umis[i].umi, umis[j].umi)
					}
				}
			}
		}
	} else {
		pairsCh := make(chan pairJob, numThreads*4)
		edgesCh := make(chan umiEdge, numThreads*4)

		var wg sync.WaitGroup
		for range numThreads {
			wg.Add(1)
			go func() {
				defer wg.Done()
				for p := range pairsCh {
					dist := umiLevenshtein(normalized[p.i], normalized[p.j])
					if dist <= umiClusterEditThreshold {
						edgesCh <- umiEdge{p.i, p.j, dist}
					}
				}
			}()
		}

		go func() {
			for i := 0; i < n; i++ {
				for j := i + 1; j < n; j++ {
					pairsCh <- pairJob{i, j}
				}
			}
			close(pairsCh)
		}()

		go func() {
			wg.Wait()
			close(edgesCh)
		}()

		for e := range edgesCh {
			edges = append(edges, e)
			if verbose {
				fmt.Printf("EDGE edit_dist=%d: %s -- %s\n", e.dist, umis[e.i].umi, umis[e.j].umi)
			}
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
	compMaxDist := make(map[int]int)
	for _, e := range edges {
		root := find(e.i)
		if e.dist > compMaxDist[root] {
			compMaxDist[root] = e.dist
		}
	}

	// Compute consensus per component via majority vote.
	compConsensus := make(map[int]string)
	for root, members := range compMembers {
		compConsensus[root] = computeConsensusUMI(members)
	}

	// Build results and populate globalConsensus map.
	results := make([]umiClusterResult, n)
	for k := range umis {
		root := find(k)
		cons := compConsensus[root]
		globalConsensus[umis[k].umi] = cons
		results[k] = umiClusterResult{
			umi:               umis[k].umi,
			consensus:         cons,
			count:             umis[k].count,
			editDist:          umiLevenshtein(normalized[k], normalizeUMISeparator(cons)),
			maxIntraClustDist: compMaxDist[root],
		}
	}
	return results
}

func countConsensus(results []umiClusterResult) int {
	seen := make(map[string]bool)
	for _, r := range results {
		seen[r.consensus] = true
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
var umiClusterBedFilename string
var umiClusterThreads int

var umiVerbose bool

func init() {
	ontUmiClusterCmd.Flags().StringVarP(&umiClusterOutput, "output", "o", "", "Output BAM file path (required)")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterTag, "umi-tag", "RX", "SAM tag containing UMI sequence")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterOrigTag, "orig-umi-tag", "OX", "SAM tag to store original UMI before correction")
	ontUmiClusterCmd.Flags().IntVar(&umiClusterOverlap, "overlap", 50, "Maximum gap (bp) between reads to group them together")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterWholeGenome, "whole-genome", false, "Process all UMIs as a single group (ignore coordinates)")
	ontUmiClusterCmd.Flags().BoolVar(&umiClusterNoStrand, "no-strand", false, "Ignore strand when grouping reads (default: group by strand)")
	ontUmiClusterCmd.Flags().BoolVarP(&umiVerbose, "verbose", "v", false, "Verbose logging")
	ontUmiClusterCmd.Flags().IntVar(&umiClusterEditThreshold, "umi-edit-distance", 3, "Maximum Levenshtein edit distance to cluster two UMIs")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterCountsFilename, "umi-counts", "", "Write UMI counts to this file")
	ontUmiClusterCmd.Flags().StringVar(&umiClusterBedFilename, "bed", "", "Write overlap regions to this BED6 file")
	ontUmiClusterCmd.Flags().IntVarP(&umiClusterThreads, "threads", "t", 1, "Threads for UMI clustering")
}
