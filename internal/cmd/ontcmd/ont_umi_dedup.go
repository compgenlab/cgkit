package ontcmd

import (
	"fmt"
	"os"
	"strconv"

	"github.com/compgen-io/cgltk/htsio"
	"github.com/spf13/cobra"
)

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
	secSupp   []*htsio.SamRecord
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

Examples:
  cgltk ont-umi-dedup -o dedup.bam --best-tag AS input.bam
  cgltk ont-umi-dedup -o dedup.bam --best-tag AS --best-tag NM- --longest input.bam
  cgltk ont-umi-dedup -o dedup.bam --best-tag AS --mark-duplicates input.bam`,
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

		return runUmiDedup(args[0], selectors)
	},
}

func runUmiDedup(inputFile string, selectors []selector) error {
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

	header.AddPGLine("ont-umi-dedup", "cgltk", "DS:UMI deduplication")

	wopts := htsio.SamWriterOptions(header).BAM().SortCoord().Threads(2)
	writer, err := htsio.NewSamWriter(umiDedupOutput, wopts)
	if err != nil {
		return err
	}

	// Active MI groups keyed by MI tag value.
	groups := make(map[string]*miGroup)

	// Track selected read names so we can handle secondary/supplementary
	// reads that arrive after their MI group has been flushed.
	selectedNames := make(map[string]bool)

	// Pending secondary/supplementary reads whose MI group hasn't been
	// flushed yet. Keyed by MI tag value.
	pendingSecSupp := make(map[string][]*htsio.SamRecord)

	currentChrom := ""

	var totalReads int64
	var totalPrimary int64
	var totalKept int64
	var totalMIGroups int64

	// flushGroup selects the best read from an MI group, writes it (and
	// its secondary/supplementary reads) to the output, and handles
	// duplicates according to --mark-duplicates.
	flushGroup := func(mi string, g *miGroup) error {
		totalMIGroups++

		if len(g.primaries) == 0 {
			// No primaries — write sec/supp through unchanged.
			for _, rec := range g.secSupp {
				if err := writer.Write(rec); err != nil {
					return err
				}
			}
			for _, rec := range pendingSecSupp[mi] {
				if err := writer.Write(rec); err != nil {
					return err
				}
			}
			delete(pendingSecSupp, mi)
			return nil
		}

		bestIdx := selectBest(g.primaries, selectors)
		bestName := g.primaries[bestIdx].ReadName
		selectedNames[bestName] = true
		totalKept++

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

		// Write secondary/supplementary: keep those belonging to the
		// selected read, discard or mark-dup the rest.
		writeSecSupp := func(recs []*htsio.SamRecord) error {
			for _, rec := range recs {
				if rec.ReadName == bestName {
					if err := writer.Write(rec); err != nil {
						return err
					}
				} else if umiDedupMarkDuplicates {
					rec.Flag |= 0x400
					if err := writer.Write(rec); err != nil {
						return err
					}
				}
			}
			return nil
		}

		if err := writeSecSupp(g.secSupp); err != nil {
			return err
		}
		if pending, ok := pendingSecSupp[mi]; ok {
			if err := writeSecSupp(pending); err != nil {
				return err
			}
			delete(pendingSecSupp, mi)
		}

		return nil
	}

	// flushExpired flushes MI groups whose maxEnd is behind curStart.
	flushExpired := func(curStart int, curChrom string) error {
		for mi, g := range groups {
			if g.chrom != curChrom || g.maxEnd <= curStart {
				if err := flushGroup(mi, g); err != nil {
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
			if err := flushGroup(mi, g); err != nil {
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
		totalReads++

		// Reads without MI tag: pass through unchanged.
		miTag, hasMI := rec.Tags[umiDedupMITag]
		if !hasMI || miTag.Value == "" {
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

		// Secondary/supplementary reads: add to their MI group if it
		// exists, otherwise check if we already know their fate.
		if rec.IsSecondary() || rec.IsSupplementary() {
			if g, ok := groups[mi]; ok {
				g.secSupp = append(g.secSupp, rec)
			} else if selectedNames[rec.ReadName] {
				// MI group already flushed, this read's primary was selected.
				if err := writer.Write(rec); err != nil {
					writer.Close()
					return err
				}
			} else if umiDedupMarkDuplicates {
				// MI group already flushed, this read's primary was NOT selected.
				rec.Flag |= 0x400
				if err := writer.Write(rec); err != nil {
					writer.Close()
					return err
				}
			}
			// else: discard
			continue
		}

		// Primary read.
		totalPrimary++

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

	fmt.Fprintf(os.Stderr, "Total reads: %d, primary: %d, MI groups: %d, kept: %d\n",
		totalReads, totalPrimary, totalMIGroups, totalKept)

	return writer.Close()
}

var umiDedupOutput string
var umiDedupBestTags []string
var umiDedupLongest bool
var umiDedupMarkDuplicates bool
var umiDedupMITag string

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
}
