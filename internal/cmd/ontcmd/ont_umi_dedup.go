package ontcmd

import (
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"

	"github.com/compgen-io/cgkit/htsio"
	"github.com/compgen-io/cgkit/htsio/bam"
	"github.com/spf13/cobra"
)

// dedupStats collects metrics during deduplication for the --stats report.
type dedupStats struct {
	totalReads      int64
	totalPrimary    int64
	totalMIGroup    int64
	totalKept       int64
	noMIReads       int64 // reads without MI tag (passed through)
	secSuppDropped  int64 // secondary/supplementary alignments removed

	// Per-MI-group sizes (number of primary reads per group).
	groupSizes []int

	// Tag values for kept and discarded primary reads, keyed by tag name.
	keptTagValues     map[string][]float64
	discardedTagValues map[string][]float64
}

func newDedupStats() *dedupStats {
	return &dedupStats{
		keptTagValues:      make(map[string][]float64),
		discardedTagValues: make(map[string][]float64),
	}
}

// recordGroup records stats for one flushed MI group.
func (s *dedupStats) recordGroup(primaries []*htsio.SamRecord, bestIdx int, tagNames []string) {
	s.totalMIGroup++
	s.groupSizes = append(s.groupSizes, len(primaries))

	for i, rec := range primaries {
		for _, tag := range tagNames {
			v, ok := numericTagValue(rec, tag)
			if !ok {
				continue
			}
			if i == bestIdx {
				s.keptTagValues[tag] = append(s.keptTagValues[tag], v)
			} else {
				s.discardedTagValues[tag] = append(s.discardedTagValues[tag], v)
			}
		}
		if i == bestIdx {
			s.totalKept++
		}
	}
}

// writeReport writes the stats report to the given path.
func (s *dedupStats) writeReport(path string) error {
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("creating stats file: %w", err)
	}
	defer f.Close()

	fmt.Fprintf(f, "=== ont-umi-dedup stats ===\n\n")

	// Summary counts.
	discarded := s.totalPrimary - s.totalKept
	dupRate := 0.0
	if s.totalPrimary > 0 {
		dupRate = float64(discarded) / float64(s.totalPrimary) * 100
	}
	fmt.Fprintf(f, "Total reads:        %d\n", s.totalReads)
	fmt.Fprintf(f, "  Primary:          %d\n", s.totalPrimary)
	fmt.Fprintf(f, "  Sec/supp dropped: %d\n", s.secSuppDropped)
	fmt.Fprintf(f, "  No MI (passthru): %d\n", s.noMIReads)
	fmt.Fprintf(f, "MI groups:          %d\n", s.totalMIGroup)
	fmt.Fprintf(f, "Reads kept:         %d\n", s.totalKept)
	fmt.Fprintf(f, "Reads discarded:    %d\n", discarded)
	fmt.Fprintf(f, "Duplication rate:   %.1f%%\n", dupRate)

	// Group size distribution.
	fmt.Fprintf(f, "\n=== reads per MI group ===\n\n")
	if len(s.groupSizes) > 0 {
		sort.Ints(s.groupSizes)
		hist := make(map[int]int)
		for _, sz := range s.groupSizes {
			hist[sz]++
		}
		// Collect and sort unique sizes.
		sizes := make([]int, 0, len(hist))
		for sz := range hist {
			sizes = append(sizes, sz)
		}
		sort.Ints(sizes)
		fmt.Fprintf(f, "  size\tcount\n")
		for _, sz := range sizes {
			fmt.Fprintf(f, "  %d\t%d\n", sz, hist[sz])
		}
		fmt.Fprintf(f, "\n")
		fmt.Fprintf(f, "  min:    %d\n", s.groupSizes[0])
		fmt.Fprintf(f, "  median: %d\n", s.groupSizes[len(s.groupSizes)/2])
		fmt.Fprintf(f, "  max:    %d\n", s.groupSizes[len(s.groupSizes)-1])
		fmt.Fprintf(f, "  mean:   %.1f\n", meanInts(s.groupSizes))
	}

	// Per-tag distributions for kept vs discarded.
	allTags := make(map[string]bool)
	for tag := range s.keptTagValues {
		allTags[tag] = true
	}
	for tag := range s.discardedTagValues {
		allTags[tag] = true
	}
	tagList := make([]string, 0, len(allTags))
	for tag := range allTags {
		tagList = append(tagList, tag)
	}
	sort.Strings(tagList)

	for _, tag := range tagList {
		fmt.Fprintf(f, "\n=== %s tag distribution ===\n\n", tag)
		kept := s.keptTagValues[tag]
		disc := s.discardedTagValues[tag]
		fmt.Fprintf(f, "  %10s  %8s  %8s  %8s\n", "", "mean", "median", "stdev")
		if len(kept) > 0 {
			fmt.Fprintf(f, "  %10s  %8.1f  %8.1f  %8.1f\n", "kept",
				meanFloats(kept), medianFloats(kept), stdevFloats(kept))
		}
		if len(disc) > 0 {
			fmt.Fprintf(f, "  %10s  %8.1f  %8.1f  %8.1f\n", "discarded",
				meanFloats(disc), medianFloats(disc), stdevFloats(disc))
		}
	}

	return nil
}

func meanInts(v []int) float64 {
	if len(v) == 0 {
		return 0
	}
	sum := 0
	for _, x := range v {
		sum += x
	}
	return float64(sum) / float64(len(v))
}

func meanFloats(v []float64) float64 {
	if len(v) == 0 {
		return 0
	}
	sum := 0.0
	for _, x := range v {
		sum += x
	}
	return sum / float64(len(v))
}

func medianFloats(v []float64) float64 {
	if len(v) == 0 {
		return 0
	}
	sorted := make([]float64, len(v))
	copy(sorted, v)
	sort.Float64s(sorted)
	return sorted[len(sorted)/2]
}

func stdevFloats(v []float64) float64 {
	if len(v) < 2 {
		return 0
	}
	m := meanFloats(v)
	sum := 0.0
	for _, x := range v {
		d := x - m
		sum += d * d
	}
	return math.Sqrt(sum / float64(len(v)))
}

// selector is a single criterion for choosing the best read from a group.
// Selectors are applied in order; each one narrows the candidates and the
// next breaks ties among the remaining reads.
type selector interface {
	// compare returns a negative value if a is better than b, positive if
	// b is better, and 0 if they are tied under this criterion.
	compare(a, b *htsio.SamRecord) int
}

// tagSelector picks the read with the best value for a numeric SAM tag.
type tagSelector struct {
	tag       string
	ascending bool // true = lower wins (e.g., NM-), false = higher wins (e.g., AS+)
}

func (s *tagSelector) compare(a, b *htsio.SamRecord) int {
	av, aOk := numericTagValue(a, s.tag)
	bv, bOk := numericTagValue(b, s.tag)
	if !aOk && !bOk {
		return 0
	}
	if !aOk {
		return 1 // a missing → b wins
	}
	if !bOk {
		return -1 // b missing → a wins
	}
	if av == bv {
		return 0
	}
	if s.ascending {
		if av < bv {
			return -1
		}
		return 1
	}
	// descending: higher wins
	if av > bv {
		return -1
	}
	return 1
}

func numericTagValue(rec *htsio.SamRecord, tag string) (float64, bool) {
	t, ok := rec.Tags[tag]
	if !ok {
		return 0, false
	}
	switch t.Type {
	case 'i':
		v, err := strconv.Atoi(t.Value)
		if err != nil {
			return 0, false
		}
		return float64(v), true
	case 'f':
		v, err := strconv.ParseFloat(t.Value, 64)
		if err != nil {
			return 0, false
		}
		return v, true
	}
	return 0, false
}

// longestSelector picks the read with the longest aligned query sequence
// (excluding soft-clipped bases).
type longestSelector struct{}

func (s *longestSelector) compare(a, b *htsio.SamRecord) int {
	al := cigarQueryAlignedLen(a.Cigar)
	bl := cigarQueryAlignedLen(b.Cigar)
	if al == bl {
		return 0
	}
	if al > bl {
		return -1
	}
	return 1
}

// cigarQueryAlignedLen returns the number of query bases consumed by a CIGAR
// string, excluding soft-clipped (S) and hard-clipped (H) bases.
// Operations M, I, =, X consume query; S consumes query but is excluded.
func cigarQueryAlignedLen(cigar string) int {
	if cigar == "*" {
		return 0
	}
	n := 0
	num := 0
	for i := 0; i < len(cigar); i++ {
		c := cigar[i]
		if c >= '0' && c <= '9' {
			num = num*10 + int(c-'0')
		} else {
			switch c {
			case 'M', 'I', '=', 'X':
				n += num
			}
			num = 0
		}
	}
	return n
}

// parseTagSelectorFlag parses a tag selector string. Supported formats:
//   - "TAG+" or "TAG" — higher value wins (+ is optional)
//   - "TAG-" — lower value wins
func parseTagSelectorFlag(s string) (*tagSelector, error) {
	if len(s) < 2 {
		return nil, fmt.Errorf("invalid selector %q: expected TAG, TAG+, or TAG-", s)
	}
	dir := s[len(s)-1]
	switch dir {
	case '+':
		return &tagSelector{tag: s[:len(s)-1], ascending: false}, nil
	case '-':
		return &tagSelector{tag: s[:len(s)-1], ascending: true}, nil
	default:
		// No suffix: default to higher wins.
		return &tagSelector{tag: s, ascending: false}, nil
	}
}

// selectBest picks the best read from a group using the ordered selectors.
// Returns the index of the best read.
func selectBest(reads []*htsio.SamRecord, selectors []selector) int {
	if len(reads) <= 1 {
		return 0
	}
	best := 0
	for i := 1; i < len(reads); i++ {
		for _, sel := range selectors {
			cmp := sel.compare(reads[i], reads[best])
			if cmp < 0 {
				best = i
				break
			}
			if cmp > 0 {
				break
			}
			// cmp == 0: tied, try next selector
		}
	}
	return best
}

// miGroup holds buffered reads for a single MI tag value.
type miGroup struct {
	primaries []*htsio.SamRecord
	maxEnd    int // max reference end (0-based exclusive) across primaries
	chrom     string
}

var ontUmiDedupCmd = &cobra.Command{
	GroupID: "ontcmd",
	Use:     "ont-umi-dedup <input.bam>",
	Short:   "Deduplicate UMI-clustered reads, keeping one representative per MI group",
	Long: `Reads a coordinate-sorted BAM file with MI tags (from ont-umi-cluster) and
selects one representative read per MI group. Selection criteria are applied
in order: each criterion narrows the candidates and the next breaks ties.

Secondary and supplementary alignments are removed from the output. In a
coordinate-sorted BAM, these alignments can appear before their primary and
cannot be reliably associated with the selected representative. Only primary
alignments are considered for selection and written to the output.

Use --threads/-t to enable parallel BGZF compression for faster output.

Examples:
  cgkit ont-umi-dedup -o dedup.bam --best-tag AS input.bam
  cgkit ont-umi-dedup -o dedup.bam --best-tag AS --best-tag NM- --longest input.bam
  cgkit ont-umi-dedup -o dedup.bam --best-tag AS --mark-duplicates input.bam
  cgkit ont-umi-dedup -o dedup.bam --best-tag AS -t 4 input.bam`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if umiDedupOutput == "" {
			return fmt.Errorf("--output is required")
		}

		// Build ordered selector chain.
		var selectors []selector
		for _, s := range umiDedupBestTags {
			ts, err := parseTagSelectorFlag(s)
			if err != nil {
				return err
			}
			selectors = append(selectors, ts)
		}
		if umiDedupLongest {
			selectors = append(selectors, &longestSelector{})
		}
		if len(selectors) == 0 {
			return fmt.Errorf("at least one selector is required (--best-tag and/or --longest)")
		}

		// Collect tag names from selectors for stats tracking.
		var statTags []string
		for _, s := range selectors {
			if ts, ok := s.(*tagSelector); ok {
				statTags = append(statTags, ts.tag)
			}
		}

		return runUmiDedup(args[0], selectors, statTags)
	},
}

func runUmiDedup(inputFile string, selectors []selector, statTags []string) error {
	reader, err := htsio.NewSamReader(inputFile)
	if err != nil {
		return err
	}
	defer reader.Close()

	header, err := reader.Header()
	if err != nil {
		return fmt.Errorf("reading header: %w", err)
	}
	if header == nil {
		return fmt.Errorf("no header found in BAM file")
	}
	if err := validateCoordinateSorted(header); err != nil {
		return err
	}

	header.AddPGLine("ont-umi-dedup", "cgkit", "DS:UMI deduplication")

	writer, err := bam.NewSortedWriter(umiDedupOutput, header, true)
	if err != nil {
		return err
	}
	if umiDedupThreads > 1 {
		writer.SetThreads(umiDedupThreads)
	}

	// Active MI groups keyed by MI tag value.
	groups := make(map[string]*miGroup)

	currentChrom := ""
	stats := newDedupStats()

	// flushGroup selects the best read from an MI group and writes it to
	// the output. Secondary/supplementary reads are dropped entirely.
	flushGroup := func(g *miGroup) error {
		if len(g.primaries) == 0 {
			return nil
		}

		bestIdx := selectBest(g.primaries, selectors)
		stats.recordGroup(g.primaries, bestIdx, statTags)

		// Write primaries: best is kept, others are discarded or marked.
		for i, rec := range g.primaries {
			if i == bestIdx {
				if err := writer.Write(rec); err != nil {
					return err
				}
			} else if umiDedupMarkDuplicates {
				rec.Flag |= 0x400
				if err := writer.Write(rec); err != nil {
					return err
				}
			}
			// else: discard (don't write)
		}

		return nil
	}

	// flushExpired flushes MI groups whose maxEnd is behind curStart.
	flushExpired := func(curStart int, curChrom string) error {
		for mi, g := range groups {
			if g.chrom != curChrom || g.maxEnd <= curStart {
				if err := flushGroup(g); err != nil {
					return err
				}
				delete(groups, mi)
			}
		}
		return nil
	}

	// flushAll flushes every remaining MI group.
	flushAll := func() error {
		for mi, g := range groups {
			if err := flushGroup(g); err != nil {
				return err
			}
			delete(groups, mi)
		}
		return nil
	}

	for rec, err := range reader.Records() {
		if err != nil {
			writer.Close()
			return err
		}
		stats.totalReads++

		// Reads without MI tag: pass through unchanged.
		miTag, hasMI := rec.Tags[umiDedupMITag]
		if !hasMI || miTag.Value == "" {
			stats.noMIReads++
			if err := writer.Write(rec); err != nil {
				writer.Close()
				return err
			}
			continue
		}
		mi := miTag.Value

		// Chromosome transition: flush all groups.
		if rec.RefName != currentChrom {
			if currentChrom != "" {
				if err := flushAll(); err != nil {
					writer.Close()
					return err
				}
			}
			currentChrom = rec.RefName
			fmt.Fprintf(os.Stderr, "Processing %s...\n", currentChrom)
		}

		// Secondary/supplementary reads are dropped entirely.
		if rec.IsSecondary() || rec.IsSupplementary() {
			stats.secSuppDropped++
			continue
		}

		// Primary read.
		stats.totalPrimary++

		readStart := rec.Pos - 1 // 0-based
		readEnd := readStart + htsio.CigarRefLen(rec.Cigar)

		// Flush expired groups before adding this read.
		if err := flushExpired(readStart, currentChrom); err != nil {
			writer.Close()
			return err
		}

		g, ok := groups[mi]
		if !ok {
			g = &miGroup{
				chrom:  rec.RefName,
				maxEnd: readEnd,
			}
			groups[mi] = g
		}
		g.primaries = append(g.primaries, rec)
		if readEnd > g.maxEnd {
			g.maxEnd = readEnd
		}
	}

	// Flush remaining groups.
	if err := flushAll(); err != nil {
		writer.Close()
		return err
	}

	fmt.Fprintf(os.Stderr, "Total reads: %d, primary: %d, sec/supp dropped: %d, MI groups: %d, kept: %d\n",
		stats.totalReads, stats.totalPrimary, stats.secSuppDropped, stats.totalMIGroup, stats.totalKept)

	if umiDedupStatsFile != "" {
		if err := stats.writeReport(umiDedupStatsFile); err != nil {
			writer.Close()
			return err
		}
	}

	return writer.Close()
}

var umiDedupOutput string
var umiDedupBestTags []string
var umiDedupLongest bool
var umiDedupMarkDuplicates bool
var umiDedupMITag string
var umiDedupStatsFile string
var umiDedupThreads int

// tagArrayValue is a pflag.Value that collects repeated --best-tag flags
// but displays as "tag" instead of "stringArray" in help text.
type tagArrayValue struct {
	values *[]string
}

func (v *tagArrayValue) String() string { return "" }
func (v *tagArrayValue) Type() string   { return "tag" }
func (v *tagArrayValue) Set(s string) error {
	*v.values = append(*v.values, s)
	return nil
}

func init() {
	ontUmiDedupCmd.Flags().StringVarP(&umiDedupOutput, "output", "o", "", "Output BAM file path (required)")
	ontUmiDedupCmd.Flags().Var(&tagArrayValue{values: &umiDedupBestTags}, "best-tag", "Tag-based selector: TAG or TAG+ (higher wins) or TAG- (lower wins); repeatable, applied in order")
	ontUmiDedupCmd.Flags().BoolVar(&umiDedupLongest, "longest", false, "Use longest query sequence as a selector (applied after --best-tag selectors)")
	ontUmiDedupCmd.Flags().BoolVar(&umiDedupMarkDuplicates, "mark-duplicates", false, "Set PCR duplicate flag (0x400) on non-selected reads instead of removing them")
	ontUmiDedupCmd.Flags().StringVar(&umiDedupMITag, "mi-tag", "MI", "SAM tag containing molecule group ID")
	ontUmiDedupCmd.Flags().StringVar(&umiDedupStatsFile, "stats", "", "Write deduplication statistics to this file")
	ontUmiDedupCmd.Flags().IntVarP(&umiDedupThreads, "threads", "t", 1, "Number of BGZF compression threads for output writing")
}
