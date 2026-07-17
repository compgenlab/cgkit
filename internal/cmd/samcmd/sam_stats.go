package samcmd

import (
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
	"text/tabwriter"

	"github.com/compgenlab/cghts/htsio"
	_ "github.com/compgenlab/cghts/htsio/bam"
	_ "github.com/compgenlab/cghts/htsio/cram"
	"github.com/spf13/cobra"
)

var (
	samStatsCramRef         string
	samStatsRgid            string
	samStatsTags            string
	samStatsCalcInsert      bool
	samStatsUnique          bool
	samStatsShowUnmappedRef bool
)

func init() {
	samStatsCmd.Flags().StringVar(&samStatsCramRef, "cram-ref", "", "Reference FASTA for CRAM files")
	samStatsCmd.Flags().StringVar(&samStatsRgid, "rgid", "", "Only count reads from this read group (RG tag)")
	samStatsCmd.Flags().StringVar(&samStatsTags, "tags", "", "Tally tag value distributions (comma-list, e.g. NH:i,RG:Z,MAPQ)")
	samStatsCmd.Flags().BoolVar(&samStatsCalcInsert, "calc-insert", false, "Calculate the median insert size of proper pairs")
	samStatsCmd.Flags().BoolVar(&samStatsUnique, "unique", false, "Only count uniquely-mapped reads")
	samStatsCmd.Flags().BoolVar(&samStatsShowUnmappedRef, "show-unmapped-ref", false, "List references even if they have zero mapped reads")
}

var samStatsCmd = &cobra.Command{
	GroupID:     "samcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "sam-stats <input>",
	Short:       "Summary statistics for a SAM/BAM/CRAM file",
	Long: `Read a SAM/BAM/CRAM file and print summary statistics: read counts,
mapping rates, Q30 percentage, depth, SAM flag breakdown, per-reference counts,
and (optionally) tag value distributions and insert-size median.`,
	Args: cobra.ExactArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		tags, err := parseTagSpecs(samStatsTags)
		if err != nil {
			return err
		}

		opts := htsio.NewSamReaderOpts()
		if samStatsCramRef != "" {
			opts.RefPath(samStatsCramRef)
		}
		reader, err := htsio.NewSamReader(args[0], opts)
		if err != nil {
			return fmt.Errorf("open %s: %w", args[0], err)
		}
		defer reader.Close()

		header, err := reader.Header()
		if err != nil {
			return fmt.Errorf("read header: %w", err)
		}

		opt := statsOpts{
			rgid:            samStatsRgid,
			tags:            tags,
			calcInsert:      samStatsCalcInsert,
			unique:          samStatsUnique,
			showUnmappedRef: samStatsShowUnmappedRef,
		}

		res, err := computeStats(reader, header, opt, cmd.ErrOrStderr())
		if err != nil {
			return err
		}

		writeStats(cmd.OutOrStdout(), res, opt)
		return nil
	},
}

// flagBits lists the SAM flag bits reported in the [Flags] section, in order.
var flagBits = []struct {
	bit  int
	name string
}{
	{0x1, "paired"},
	{0x2, "proper pair"},
	{0x4, "unmapped"},
	{0x8, "mate unmapped"},
	{0x10, "reverse strand"},
	{0x20, "mate reverse strand"},
	{0x40, "first in pair"},
	{0x80, "second in pair"},
	{0x100, "secondary"},
	{0x200, "QC fail"},
	{0x400, "duplicate"},
	{0x800, "supplementary"},
}

// tagSpec describes a tag to tally. name is the two-char tag (or "MAPQ");
// numeric controls whether the distribution is sorted numerically.
type tagSpec struct {
	name    string
	numeric bool
}

// parseTagSpecs parses a "--tags" value like "NH:i,RG:Z,MAPQ".
func parseTagSpecs(s string) ([]tagSpec, error) {
	if strings.TrimSpace(s) == "" {
		return nil, nil
	}
	var specs []tagSpec
	for _, raw := range strings.Split(s, ",") {
		k := strings.TrimSpace(raw)
		if k == "" {
			continue
		}
		if strings.EqualFold(k, "MAPQ") {
			specs = append(specs, tagSpec{name: "MAPQ", numeric: true})
			continue
		}
		idx := strings.Index(k, ":")
		if idx < 0 {
			return nil, fmt.Errorf("invalid tag %q: specify a type (e.g. NH:i, RG:Z) or the special MAPQ", k)
		}
		name := k[:idx]
		typ := strings.ToUpper(k[idx+1:])
		switch typ {
		case "Z":
			specs = append(specs, tagSpec{name: name, numeric: false})
		case "I":
			specs = append(specs, tagSpec{name: name, numeric: true})
		default:
			return nil, fmt.Errorf("invalid tag %q: only Z (string) and I (integer) types are supported", k)
		}
	}
	return specs, nil
}

type statsOpts struct {
	rgid            string
	tags            []tagSpec
	calcInsert      bool
	unique          bool
	showUnmappedRef bool
}

type statsResult struct {
	total      int
	mapped     int
	unmapped   int
	multiple   int
	totalBases int64
	q30Bases   int64
	refLength  int64
	hasGaps    bool

	flagCounts map[int]int

	refCounts map[string]int
	refOrder  []string

	tagCounts  map[string]map[string]int
	tagMissing map[string]int

	insertSizes      []int
	paired           bool
	insertDisabled   bool // calc-insert disabled because gapped alignments were seen
	calcInsertActive bool // calc-insert requested and still active at end
}

// isUniquelyMapped ports ReadUtils.isReadUniquelyMapped: prefer the NH tag,
// then IH, falling back to a non-zero MAPQ.
func isUniquelyMapped(rec *htsio.SamRecord) bool {
	if t, ok := rec.Tags["NH"]; ok {
		if v, ok := t.IntValue(); ok {
			return v == 1
		}
	}
	if t, ok := rec.Tags["IH"]; ok {
		if v, ok := t.IntValue(); ok {
			return v == 1
		}
	}
	return rec.MapQ != 0
}

// countBases walks a CIGAR and returns the number of query bases consumed by
// M/=/X/I operations, how many of those have phred quality >= 30, and whether
// any gapped (N) operation was present.
func countBases(cigar, qual string) (bases, q30 int, hasGaps bool) {
	if cigar == "*" || cigar == "" {
		return 0, 0, false
	}
	hasQual := qual != "*" && qual != ""
	qpos := 0
	num := 0
	for i := 0; i < len(cigar); i++ {
		c := cigar[i]
		if c >= '0' && c <= '9' {
			num = num*10 + int(c-'0')
			continue
		}
		switch c {
		case 'M', '=', 'X', 'I':
			bases += num
			if hasQual {
				for j := 0; j < num; j++ {
					// Phred+33 encoding: quality 30 is ASCII 63 ('?').
					if qpos+j < len(qual) && qual[qpos+j] >= 63 {
						q30++
					}
				}
			}
			qpos += num
		case 'S':
			qpos += num
		case 'N':
			hasGaps = true
		}
		num = 0
	}
	return bases, q30, hasGaps
}

func computeStats(reader htsio.SamReader, header *htsio.SamHeader, opt statsOpts, warnw io.Writer) (*statsResult, error) {
	res := &statsResult{
		flagCounts:       make(map[int]int),
		refCounts:        make(map[string]int),
		tagCounts:        make(map[string]map[string]int),
		tagMissing:       make(map[string]int),
		calcInsertActive: opt.calcInsert,
	}
	for _, spec := range opt.tags {
		res.tagCounts[spec.name] = make(map[string]int)
	}
	for _, ref := range header.References() {
		res.refOrder = append(res.refOrder, ref.Name)
		res.refCounts[ref.Name] = 0
		res.refLength += int64(ref.Length)
	}

	calcInsert := opt.calcInsert

	for rec, err := range reader.Records() {
		if err != nil {
			return nil, fmt.Errorf("read record: %w", err)
		}

		// --rgid: restrict to a single read group.
		if opt.rgid != "" {
			rg, ok := rec.Tags["RG"]
			if !ok || rg.Value != opt.rgid {
				continue
			}
		}

		// Base / Q30 accounting: non-duplicate, mapped, and unique-passing reads.
		if !rec.IsDuplicate() && !rec.IsUnmapped() && (isUniquelyMapped(rec) || !opt.unique) {
			bases, q30, hasGaps := countBases(rec.Cigar, rec.Qual)
			res.totalBases += int64(bases)
			res.q30Bases += int64(q30)
			if hasGaps {
				res.hasGaps = true
				if calcInsert {
					calcInsert = false
					res.insertDisabled = true
					res.calcInsertActive = false
					fmt.Fprintln(warnw, "Warning: gapped alignments found - not calculating insert size.")
				}
			}
		}

		// Profile only the first read of a pair from here on.
		if rec.Flag&0x1 != 0 && rec.Flag&0x80 != 0 {
			continue
		}

		res.total++

		for _, fb := range flagBits {
			if rec.Flag&fb.bit != 0 {
				res.flagCounts[fb.bit]++
			}
		}

		// Skip duplicates, QC failures, and supplementary alignments.
		if rec.IsDuplicate() || rec.Flag&0x200 != 0 || rec.IsSupplementary() {
			continue
		}

		if rec.IsUnmapped() {
			res.unmapped++
			continue
		}
		res.mapped++
		if !isUniquelyMapped(rec) {
			res.multiple++
			if opt.unique {
				continue
			}
		}

		// Tag value distributions.
		for _, spec := range opt.tags {
			if spec.name == "MAPQ" {
				res.tagCounts["MAPQ"][strconv.Itoa(rec.MapQ)]++
				continue
			}
			if t, ok := rec.Tags[spec.name]; ok {
				res.tagCounts[spec.name][t.Value]++
			} else {
				res.tagMissing[spec.name]++
			}
		}

		// Insert size: proper pairs whose mate is on the same reference.
		if calcInsert && rec.Flag&0x1 != 0 {
			res.paired = true
			if rec.Flag&0x2 != 0 && (rec.RefNext == "=" || rec.RefNext == rec.RefName) {
				size := rec.InsertLen
				if size < 0 {
					size = -size
				}
				res.insertSizes = append(res.insertSizes, size)
			}
		}

		res.refCounts[rec.RefName]++
	}

	return res, nil
}

func median(vals []int) int {
	if len(vals) == 0 {
		return 0
	}
	s := append([]int(nil), vals...)
	sort.Ints(s)
	n := len(s)
	if n%2 == 1 {
		return s[n/2]
	}
	return (s[n/2-1] + s[n/2]) / 2
}

func writeStats(w io.Writer, res *statsResult, opt statsOpts) {
	tw := tabwriter.NewWriter(w, 0, 0, 2, ' ', 0)
	fmt.Fprintf(tw, "Total reads:\t%d\n", res.total)
	fmt.Fprintf(tw, "Mapped reads:\t%d\n", res.mapped)
	fmt.Fprintf(tw, "Unmapped reads:\t%d\n", res.unmapped)
	fmt.Fprintf(tw, "Multiple-mapped reads:\t%d\n", res.multiple)
	fmt.Fprintf(tw, "Uniquely-mapped reads:\t%d\n", res.mapped-res.multiple)
	if res.totalBases > 0 {
		fmt.Fprintf(tw, "Q30 pct:\t%.2f%%\n", 100.0*float64(res.q30Bases)/float64(res.totalBases))
	} else {
		fmt.Fprintf(tw, "Q30 pct:\t0.00%%\n")
	}
	fmt.Fprintf(tw, "Total bases:\t%d\n", res.totalBases)
	fmt.Fprintf(tw, "Reference length:\t%d\n", res.refLength)
	if !res.hasGaps && res.refLength > 0 {
		fmt.Fprintf(tw, "Total depth:\t%.2fX\n", float64(res.totalBases)/float64(res.refLength))
	}
	if res.calcInsertActive && res.paired {
		fmt.Fprintf(tw, "Median insert size:\t%d\n", median(res.insertSizes))
	}
	tw.Flush()

	// [Flags]
	fmt.Fprintln(w)
	fmt.Fprintln(w, "[Flags]")
	ftw := tabwriter.NewWriter(w, 0, 0, 2, ' ', 0)
	for _, fb := range flagBits {
		if c := res.flagCounts[fb.bit]; c > 0 {
			fmt.Fprintf(ftw, "0x%x\t%s\t%d\n", fb.bit, fb.name, c)
		}
	}
	ftw.Flush()

	// Tag distributions.
	for _, spec := range opt.tags {
		fmt.Fprintln(w)
		fmt.Fprintf(w, "[%s]\n", spec.name)
		ttw := tabwriter.NewWriter(w, 0, 0, 2, ' ', 0)
		counts := res.tagCounts[spec.name]
		keys := make([]string, 0, len(counts))
		for k := range counts {
			keys = append(keys, k)
		}
		if spec.numeric {
			sort.Slice(keys, func(i, j int) bool {
				ai, _ := strconv.Atoi(keys[i])
				aj, _ := strconv.Atoi(keys[j])
				return ai < aj
			})
		} else {
			sort.Strings(keys)
		}
		for _, k := range keys {
			fmt.Fprintf(ttw, "%s\t%d\n", k, counts[k])
		}
		if m := res.tagMissing[spec.name]; m > 0 {
			fmt.Fprintf(ttw, "missing\t%d\n", m)
		}
		ttw.Flush()
	}

	// [References]
	fmt.Fprintln(w)
	fmt.Fprintln(w, "[References]")
	rtw := tabwriter.NewWriter(w, 0, 0, 2, ' ', 0)
	for _, ref := range res.refOrder {
		c := res.refCounts[ref]
		if c == 0 && !opt.showUnmappedRef {
			continue
		}
		fmt.Fprintf(rtw, "%s\t%d\n", ref, c)
	}
	rtw.Flush()
}
