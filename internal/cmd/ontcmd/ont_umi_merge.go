package ontcmd

import (
	"fmt"
	"io"
	"os"
	"sort"
	"strings"

	"github.com/compgen-io/cgltk/align"
	"github.com/compgen-io/cgltk/htsio"
	"github.com/compgen-io/cgltk/seqio"
	"github.com/spf13/cobra"
)

var ontUmiMergeCmd = &cobra.Command{
	GroupID: "ontcmd",
	Use:     "ont-umi-merge <input.bam>",
	Short:   "Collapse similar UMIs in a coordinate-sorted BAM file",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if umiMergeOutput == "" {
			return fmt.Errorf("--output is required")
		}

		inputFile := args[0]

		// Open counts writer if requested
		var countsWriter io.Writer
		var closeCounts func() error
		if umiMergeCountsFilename != "" {
			var err error
			countsWriter, closeCounts, err = openWriter(umiMergeCountsFilename)
			if err != nil {
				return fmt.Errorf("opening umi-counts: %w", err)
			}
			defer func() {
				if err := closeCounts(); err != nil {
					fmt.Fprintf(os.Stderr, "error closing umi-counts: %v\n", err)
				}
			}()
		}

		if umiMergeWholeGenome {
			return umiMergeWholeGenomeMode(inputFile, countsWriter)
		}
		return umiMergeOverlapMode(inputFile, countsWriter)
	},
}

// umiMergeOverlapMode processes each overlap group independently:
// buffer the group's reads, cluster UMIs, update tags, write out, repeat.
func umiMergeOverlapMode(inputFile string, countsWriter io.Writer) error {
	reader := htsio.NewSamReader(inputFile)

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

	header.AddLine(fmt.Sprintf("@PG\tID:ont-umi-merge\tPN:cgltk\tCL:ont-umi-merge\tDS:UMI collapsing; canonical UMI written to %s, original preserved in %s", umiMergeTag, umiMergeOrigTag))

	writer := htsio.NewSamWriter(umiMergeOutput, header).Format(htsio.FormatBAM)

	// Current overlap group state
	var groupRecs []*htsio.SamRecord
	var groupRName string
	var groupMinStart int // 0-based
	var groupMaxEnd int   // 0-based, exclusive

	totalReads := 0
	totalChanged := 0

	flushGroup := func() error {
		if len(groupRecs) == 0 {
			return nil
		}

		// Collect UMI counts for this group
		umiCounts := make(map[string]int)
		for _, rec := range groupRecs {
			if umi := getUMI(rec); umi != "" {
				umiCounts[umi]++
			}
		}

		// Cluster UMIs
		canonical := make(map[string]string)
		results := clusterUMIs(umiCounts, canonical, umiVerbose)

		// Write counts
		writeGroupCounts(countsWriter, results, groupRName, groupMinStart, groupMaxEnd, false)

		canonicalCount := countCanonical(results)
		fmt.Fprintf(os.Stderr, "%s:%d-%d: %d reads, %d unique UMIs -> %d canonical\n",
			groupRName, groupMinStart, groupMaxEnd, len(groupRecs), len(umiCounts), canonicalCount)

		// Update tags and write reads
		for _, rec := range groupRecs {
			tag, hasTag := rec.Tags[umiMergeTag]
			if hasTag {
				umi := tag.Value
				if can, ok := canonical[umi]; ok && can != umi {
					rec.Tags[umiMergeOrigTag] = htsio.SamTag{Type: 'Z', Value: umi}
					rec.Tags[umiMergeTag] = htsio.SamTag{Type: 'Z', Value: can}
					totalChanged++
				}
			}
			totalReads++
			if err := writer.Write(rec); err != nil {
				return err
			}
		}

		groupRecs = nil
		return nil
	}

	addToGroup := func(rec *htsio.SamRecord) error {
		readStart := rec.Pos - 1 // convert to 0-based
		readEnd := readStart + htsio.CigarRefLen(rec.Cigar)

		if rec.IsUnmapped() || rec.Cigar == "*" {
			// Unmapped reads: flush current group, write through immediately
			if err := flushGroup(); err != nil {
				return err
			}
			totalReads++
			return writer.Write(rec)
		}

		if len(groupRecs) == 0 {
			// Start new group
			groupRName = rec.RName
			groupMinStart = readStart
			groupMaxEnd = readEnd
			groupRecs = append(groupRecs, rec)
			return nil
		}

		// Check overlap with current group
		overlap := min(groupMaxEnd, readEnd) - readStart
		if rec.RName != groupRName || overlap < umiMergeOverlap {
			// Flush current group, start new one
			if err := flushGroup(); err != nil {
				return err
			}
			groupRName = rec.RName
			groupMinStart = readStart
			groupMaxEnd = readEnd
			groupRecs = append(groupRecs, rec)
			return nil
		}

		// Extends current group
		if readEnd > groupMaxEnd {
			groupMaxEnd = readEnd
		}
		groupRecs = append(groupRecs, rec)
		return nil
	}

	if firstRec != nil {
		if err := addToGroup(firstRec); err != nil {
			return err
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
		if err := addToGroup(rec); err != nil {
			return err
		}
	}
	reader.Close()

	// Flush final group
	if err := flushGroup(); err != nil {
		return err
	}

	fmt.Fprintf(os.Stderr, "Total reads: %d, UMIs corrected: %d\n", totalReads, totalChanged)
	return writer.Close()
}

// umiMergeWholeGenomeMode uses two passes over the entire file:
// pass 1 collects all UMI counts, pass 2 rewrites with canonical UMIs.
func umiMergeWholeGenomeMode(inputFile string, countsWriter io.Writer) error {
	// Pass 1: collect all UMIs
	reader := htsio.NewSamReader(inputFile)
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
	globalCanonical := make(map[string]string)
	results := clusterUMIs(umiCounts, globalCanonical, umiVerbose)
	writeGroupCounts(countsWriter, results, "", 0, 0, true)
	canonicalCount := countCanonical(results)
	fmt.Fprintf(os.Stderr, "whole-genome: %d unique UMIs -> %d canonical\n", len(umiCounts), canonicalCount)

	// Pass 2: rewrite BAM
	header.AddLine(fmt.Sprintf("@PG\tID:ont-umi-merge\tPN:cgltk\tCL:ont-umi-merge\tDS:UMI collapsing; canonical UMI written to %s, original preserved in %s", umiMergeTag, umiMergeOrigTag))

	reader2 := htsio.NewSamReader(inputFile)
	writer := htsio.NewSamWriter(umiMergeOutput, header).Format(htsio.FormatBAM)

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

		tag, hasTag := rec.Tags[umiMergeTag]
		if hasTag {
			umi := tag.Value
			if canonical, ok := globalCanonical[umi]; ok && canonical != umi {
				rec.Tags[umiMergeOrigTag] = htsio.SamTag{Type: 'Z', Value: umi}
				rec.Tags[umiMergeTag] = htsio.SamTag{Type: 'Z', Value: canonical}
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

func writeGroupCounts(w io.Writer, results []umiClusterResult, rname string, start0, end0 int, wholeGenome bool) {
	if w == nil || len(results) == 0 {
		return
	}

	// Sum counts per canonical UMI
	canonicalTotals := make(map[string]int)
	for _, r := range results {
		canonicalTotals[r.canonical] += r.count
	}

	for _, r := range results {
		if !wholeGenome {
			fmt.Fprintf(w, "%s\t%d\t%d\t", rname, start0, end0)
		}
		fmt.Fprintf(w, "%s\t%s\t%d\t%d\t%d\n", r.umi, r.canonical, r.count, canonicalTotals[r.canonical], r.matchScore)
	}
}

func getUMI(rec *htsio.SamRecord) string {
	tag, ok := rec.Tags[umiMergeTag]
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

// detectSeparator returns the UMI separator: "-" or "TT".
func detectSeparator(umi string) string {
	if strings.Contains(umi, "-") {
		return "-"
	}
	return "TT"
}

// countNonSepBases returns the number of non-separator bases in a UMI.
func countNonSepBases(umi string, sep string) int {
	if sep == "-" {
		return len(umi) - strings.Count(umi, "-")
	}
	// For TT separator, count non-T bases.
	// UMI pattern is VVVVTTVVVVTTVVVVTTVVVV where V is non-T.
	// Since UMI code bases are A,C,G (the V positions from the VVVV pattern),
	// all T's in the UMI string are separators.
	count := 0
	for i := 0; i < len(umi); i++ {
		if umi[i] != 'T' {
			count++
		}
	}
	return count
}

// umiClusterResult holds the clustering result for one UMI.
type umiClusterResult struct {
	umi        string
	canonical  string
	count      int
	matchScore int // non-separator matches vs canonical
}

type umiCount struct {
	umi   string
	count int
}

func clusterUMIs(umiCounts map[string]int, globalCanonical map[string]string, verbose bool) []umiClusterResult {
	if len(umiCounts) <= 1 {
		// Single UMI or empty, nothing to cluster
		var results []umiClusterResult
		for umi, count := range umiCounts {
			globalCanonical[umi] = umi
			sep := detectSeparator(umi)
			results = append(results, umiClusterResult{umi, umi, count, countNonSepBases(umi, sep)})
		}
		return results
	}

	// Sort by count descending
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

	// Detect separator from first UMI
	sep := detectSeparator(umis[0].umi)

	// clusterOf[i] = index of the canonical UMI for cluster containing i
	clusterOf := make([]int, len(umis))
	matchScores := make([]int, len(umis)) // non-sep matches vs canonical
	for i := range clusterOf {
		clusterOf[i] = -1
		matchScores[i] = countNonSepBases(umis[i].umi, sep) // self-match = all non-sep bases
	}

	opts := align.OntAlignmentDefaults() //.ClippingDisable()
	aligner := align.NewGlobalAligner(opts)

	for i := 0; i < len(umis); i++ {
		if clusterOf[i] > -1 {
			continue // already merged
		}

		// We are sorted desc, so if we are looking at this UMI in the
		// first loop, this is an anchor UMI/seq

		clusterOf[i] = i
		seqA := seqio.NewStringSeq(umis[i].umi, "a").FullSeq()

		for j := i + 1; j < len(umis); j++ {
			if clusterOf[j] > -1 {
				continue // already merged
			}

			seqB := seqio.NewStringSeq(umis[j].umi, "b").FullSeq()
			aln := aligner.Align(seqA, seqB)

			nonSepMatches := countNonSepAlignedMatches(aln, sep)
			if nonSepMatches >= umiMergeMatchThreshold {
				clusterOf[j] = i
				matchScores[j] = nonSepMatches

				if verbose {
					fmt.Printf("MATCH! %d\n%s\n===\n", nonSepMatches, aln.String())
				}
			}
		}
	}

	// Build canonical mapping and results
	results := make([]umiClusterResult, len(umis))
	for k := range umis {
		canonical := umis[clusterOf[k]].umi
		globalCanonical[umis[k].umi] = canonical
		results[k] = umiClusterResult{umis[k].umi, canonical, umis[k].count, matchScores[k]}
	}
	return results
}

// countNonSepAlignedMatches counts matching bases at non-separator positions
// in the alignment.
func countNonSepAlignedMatches(aln *align.PairwiseAlignment, sep string) int {
	qAligned := aln.QueryAlignedStr()
	tAligned := aln.TargetAlignedStr()

	matches := 0
	for i := 0; i < len(qAligned) && i < len(tAligned); i++ {
		q := qAligned[i]
		t := tAligned[i]
		if q == '-' || t == '-' {
			continue
		}
		// Skip separator positions
		if isSeparatorChar(q, sep) {
			continue
		}
		if q == t {
			matches++
		}
	}
	return matches
}

func countCanonical(results []umiClusterResult) int {
	seen := make(map[string]bool)
	for _, r := range results {
		seen[r.canonical] = true
	}
	return len(seen)
}

func isSeparatorChar(c byte, sep string) bool {
	if sep == "-" {
		return c == '-'
	}
	// TT separator: T's are separator chars
	return c == 'T'
}

var umiMergeOutput string
var umiMergeTag string
var umiMergeOrigTag string
var umiMergeOverlap int
var umiMergeWholeGenome bool
var umiMergeMatchThreshold int
var umiMergeCountsFilename string

var umiVerbose bool

func init() {
	ontUmiMergeCmd.Flags().StringVarP(&umiMergeOutput, "output", "o", "", "Output BAM file path (required)")
	ontUmiMergeCmd.Flags().StringVar(&umiMergeTag, "umi-tag", "RX", "SAM tag containing UMI sequence")
	ontUmiMergeCmd.Flags().StringVar(&umiMergeOrigTag, "orig-umi-tag", "OX", "SAM tag to store original UMI before correction")
	ontUmiMergeCmd.Flags().IntVar(&umiMergeOverlap, "overlap", 50, "Minimum overlap (bp) to group reads together")
	ontUmiMergeCmd.Flags().BoolVar(&umiMergeWholeGenome, "whole-genome", false, "Process all UMIs as a single group (ignore coordinates)")
	ontUmiMergeCmd.Flags().BoolVarP(&umiVerbose, "verbose", "v", false, "Verbose logging")
	ontUmiMergeCmd.Flags().IntVar(&umiMergeMatchThreshold, "umi-match-threshold", 13, "Minimum matching non-separator bases to cluster two UMIs")
	ontUmiMergeCmd.Flags().StringVar(&umiMergeCountsFilename, "umi-counts", "", "Write UMI counts to this file")
}
