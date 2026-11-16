package samcmd

import (
	"fmt"
	"io"
	"math"
	"math/big"
	"sort"
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/bed"
	"github.com/compgenlab/cghts/gtf"
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
	samStatsGtf             string
	samStatsLenient         bool
	samStatsSilent          bool
)

func init() {
	samStatsCmd.Flags().StringVar(&samStatsCramRef, "cram-ref", "", "Reference FASTA for CRAM files (path, http(s):// URL, or s3://)")
	samStatsCmd.Flags().StringVar(&samStatsRgid, "rgid", "", "Only count reads from this read group (RG tag)")
	samStatsCmd.Flags().StringVar(&samStatsTags, "tags", "", "Tally tag value distributions (comma-list, e.g. NH:i,RG:Z,MAPQ)")
	samStatsCmd.Flags().BoolVar(&samStatsCalcInsert, "calc-insert", false, "Calculate the median insert size of proper pairs")
	samStatsCmd.Flags().BoolVar(&samStatsUnique, "unique", false, "Only count uniquely-mapped reads")
	samStatsCmd.Flags().BoolVar(&samStatsShowUnmappedRef, "show-unmapped-ref", false, "List references even if they have zero mapped reads")
	samStatsCmd.Flags().StringVar(&samStatsGtf, "gtf", "", "GTF annotation: add gene-region counts and the library orientation (whole file loaded into memory)")
	// Accepted so pipelines written for ngsutilsj bam-stats run unchanged.
	// htsjdk validates records strictly unless told otherwise; cghts does not,
	// so both are no-ops here.
	samStatsCmd.Flags().BoolVar(&samStatsLenient, "lenient", false, "Accepted for ngsutilsj compatibility (no effect: validation is already lenient)")
	samStatsCmd.Flags().BoolVar(&samStatsSilent, "silent", false, "Accepted for ngsutilsj compatibility (no effect: validation is already lenient)")
}

var samStatsCmd = &cobra.Command{
	GroupID:     "samcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "sam-stats <input>",
	Short:       "Summary statistics for a SAM/BAM/CRAM file",
	Long: `Read a SAM/BAM/CRAM file and print summary statistics: read counts,
mapping rates, Q30 percentage, depth, SAM flag breakdown, per-reference counts,
and (optionally) tag value distributions, insert-size median and gene-region
counts.

The report is a drop-in replacement for ngsutilsj bam-stats: the same sections,
labels, ordering and number formatting, so a pipeline that parses one parses
the other. Only the first read of a pair is profiled.

--gtf classifies each counted read (unpaired, or the first read of a proper
pair) by the gene region at its 5' end, and from the sense/anti-sense balance
reports the library orientation (FR, RF or unstranded), exonic, intronic and
junction counts and their ratios. The GTF is loaded whole into memory, as
ngsutilsj does: a few GB for GENCODE.

Gene-region counts can differ slightly from ngsutilsj's. Where a position's
transcripts or overlapping genes disagree, ngsutilsj's answer depends on the
order it visits them (Java hash order, partly JVM identity hashes); cghts
classifies a position the same way whatever that order. Everything else in the
report matches ngsutilsj exactly.`,
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
		reader, err := htsio.OpenSamReader(cmd.Context(), args[0], opts)
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
		if samStatsGtf != "" {
			genes, closeGenes, err := openGeneModel(cmd.ErrOrStderr(), samStatsGtf)
			if err != nil {
				return err
			}
			defer closeGenes()
			opt.genes = genes
		}

		res, err := computeStats(reader, header, opt, cmd.ErrOrStderr())
		if err != nil {
			return err
		}

		writeStats(cmd.OutOrStdout(), res, opt)
		return nil
	},
}

// flagBits lists the SAM flag bits reported in the [Flags] section, in order,
// with ngsutilsj's labels -- including its "complimented" -- since pipelines
// match on them.
var flagBits = []struct {
	bit  int
	name string
}{
	{0x1, "Multiple fragments"},
	{0x2, "All fragments aligned"},
	{0x4, "Unmapped"},
	{0x8, "Next unmapped"},
	{0x10, "Reverse complimented"},
	{0x20, "Next reverse complimented"},
	{0x40, "First fragment"},
	{0x80, "Last fragment"},
	{0x100, "Secondary alignment"},
	{0x200, "QC Fail"},
	{0x400, "PCR/Optical duplicate"},
	{0x800, "Supplementary"},
}

// tagSpec describes a tag to tally. name is the two-char tag (or "MAPQ");
// numeric controls whether it is tallied as integers or strings.
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

// openGeneModel loads the --gtf annotation whole into memory, as ngsutilsj
// does. cghts can also classify through a tabix-indexed GTF without loading it,
// but that path costs tens of milliseconds a read today -- its region reads scan
// forward from the first index chunk -- which on a BAM of millions of reads is
// hours. Revisit when that is fixed; the classification is identical.
func openGeneModel(warnw io.Writer, path string) (geneModel, func() error, error) {
	fmt.Fprint(warnw, "Loading GTF annotation: "+path+"...")
	src, err := gtf.NewAnnotationSource(path, nil)
	if err != nil {
		fmt.Fprintln(warnw)
		return nil, nil, err
	}
	fmt.Fprintln(warnw, " done")
	return src, func() error { return nil }, nil
}

// geneModel is what the --gtf counter needs from a gene annotation.
type geneModel interface {
	FindGenicRegionForPos(ref string, pos int, strand bed.Strand, geneID string) gtf.GenicRegion
}

type statsOpts struct {
	rgid            string
	tags            []tagSpec
	calcInsert      bool
	unique          bool
	showUnmappedRef bool
	genes           geneModel // nil without --gtf
}

// intTally ports compgen-common TallyCounts: integer values, negative ones
// ignored outright (neither counted nor missing), written as every integer
// from min to max -- zeros included -- after the missing count. With nothing
// counted min and max stay -1, so an all-missing tag writes "-1\t0".
type intTally struct {
	counts   map[int]int64
	min, max int
	total    int64
	missing  int64
}

func newIntTally() *intTally { return &intTally{counts: map[int]int64{}, min: -1, max: -1} }

func (t *intTally) incr(k int) {
	if k < 0 {
		return
	}
	t.total++
	if _, ok := t.counts[k]; !ok {
		if t.min == -1 || k < t.min {
			t.min = k
		}
		if t.max == -1 || k > t.max {
			t.max = k
		}
	}
	t.counts[k]++
}

// median ports TallyCounts.getQuantile(0.5), quirk included: it is the first
// value whose cumulative count exceeds half the total, and when a value's
// count lands exactly on half the loop stops and the maximum is returned.
func (t *intTally) median() int {
	thres := 0.5 * float64(t.total)
	var count int64
	for i := t.min; i <= t.max && float64(count) < thres; i++ {
		count += t.counts[i]
		if float64(count) > thres {
			return i
		}
	}
	return t.max
}

func (t *intTally) write(w io.Writer) {
	if t.missing > 0 {
		fmt.Fprintf(w, "missing\t%d\n", t.missing)
	}
	for i := t.min; i <= t.max; i++ {
		fmt.Fprintf(w, "%d\t%d\n", i, t.counts[i])
	}
}

// strTally ports compgen-common TallyValues<String>: values in sorted order
// after the missing count.
type strTally struct {
	counts  map[string]int64
	missing int64
}

func (t *strTally) write(w io.Writer) {
	if t.missing > 0 {
		fmt.Fprintf(w, "missing\t%d\n", t.missing)
	}
	keys := make([]string, 0, len(t.counts))
	for k := range t.counts {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, k := range keys {
		fmt.Fprintf(w, "%s\t%d\n", k, t.counts[k])
	}
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
	refOrder  []string // header order: the insertion order of ngsutilsj's HashMap

	strTags map[string]*strTally
	intTags map[string]*intTally
	// Each kind in the order its keys were first given, which with
	// javaHashOrder reproduces ngsutilsj's HashMap iteration.
	strOrder, intOrder []string

	inserts *intTally
	paired  bool

	calcInsertActive bool // calc-insert requested and still active at end

	regions map[gtf.GenicRegion]int64 // nil without --gtf
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

// readGenicRegion ports GTFAnnotationSource.findGenicRegion(read, FR): the
// region at the read's 5' end, on the strand an FR library puts it, made a
// junction when the read is spliced and lands in a gene.
//
// The position is ngsutilsj's, quirk included: a forward read's first aligned
// base, but for a reverse read getAlignmentEnd() -- 1-based inclusive -- used as
// a 0-based coordinate, i.e. the base just past the read. Correcting it would
// change which reads near an exon edge count as exonic.
func readGenicRegion(genes geneModel, rec *htsio.SamRecord) gtf.GenicRegion {
	pos, strand := rec.AlignStart(), bed.StrandPlus
	if rec.IsReverse() {
		pos, strand = rec.AlignEnd(), bed.StrandMinus
	}
	reg := genes.FindGenicRegionForPos(rec.RefName, pos, strand, "")
	if strings.ContainsRune(rec.Cigar, 'N') {
		return gtf.Junctionize(reg)
	}
	return reg
}

func computeStats(reader htsio.SamReader, header *htsio.SamHeader, opt statsOpts, warnw io.Writer) (*statsResult, error) {
	res := &statsResult{
		flagCounts:       make(map[int]int),
		refCounts:        make(map[string]int),
		strTags:          map[string]*strTally{},
		intTags:          map[string]*intTally{},
		inserts:          newIntTally(),
		calcInsertActive: opt.calcInsert,
	}
	for _, spec := range opt.tags {
		if spec.numeric {
			if _, ok := res.intTags[spec.name]; !ok {
				res.intOrder = append(res.intOrder, spec.name)
			}
			res.intTags[spec.name] = newIntTally()
		} else {
			if _, ok := res.strTags[spec.name]; !ok {
				res.strOrder = append(res.strOrder, spec.name)
			}
			res.strTags[spec.name] = &strTally{counts: map[string]int64{}}
		}
	}
	if opt.genes != nil {
		res.regions = map[gtf.GenicRegion]int64{}
	}
	for _, ref := range header.References() {
		res.refOrder = append(res.refOrder, ref.Name)
		res.refCounts[ref.Name] = 0
		res.refLength += int64(ref.Length)
	}
	// --rgid matches through the header, as htsjdk's getReadGroup() does: a
	// read whose RG is not declared by an @RG line has no read group at all.
	declaredRG := map[string]bool{}
	for _, id := range header.ReadGroups() {
		declaredRG[id] = true
	}

	calcInsert := opt.calcInsert

	for rec, err := range reader.Records() {
		if err != nil {
			return nil, fmt.Errorf("read record: %w", err)
		}

		if opt.rgid != "" {
			rg, ok := rec.Tags["RG"]
			if !ok || !declaredRG[rg.Value] || rg.Value != opt.rgid {
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
					res.calcInsertActive = false
					fmt.Fprintln(warnw, "Warning: Gapped alignments found - not calculating insert-size.")
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

		for _, name := range res.strOrder {
			t, ok := rec.Tags[name]
			if !ok {
				res.strTags[name].missing++
				continue
			}
			v, ok := t.StringValue()
			if !ok {
				return nil, fmt.Errorf("read %s: tag %s is type %c, not a string (Z)", rec.ReadName, name, t.Type)
			}
			res.strTags[name].counts[v]++
		}
		for _, name := range res.intOrder {
			if name == "MAPQ" {
				res.intTags[name].incr(rec.MapQ)
				continue
			}
			t, ok := rec.Tags[name]
			if !ok {
				res.intTags[name].missing++
				continue
			}
			v, ok := t.IntValue()
			if !ok {
				return nil, fmt.Errorf("read %s: tag %s is type %c, not an integer", rec.ReadName, name, t.Type)
			}
			res.intTags[name].incr(v)
		}

		// Insert size: proper pairs whose mate is on the same reference.
		if calcInsert && rec.Flag&0x1 != 0 {
			res.paired = true
			if rec.Flag&0x2 != 0 && (rec.RefNext == "=" || rec.RefNext == rec.RefName) {
				size := rec.InsertLen
				if size < 0 {
					size = -size
				}
				res.inserts.incr(size)
			}
		}

		res.refCounts[rec.RefName]++

		// Gene regions: unpaired reads, and the first read of a proper pair.
		// (Duplicates, QC failures, supplementary and -- under --unique --
		// multi-mapped reads were already skipped above.)
		if res.regions != nil && (rec.Flag&0x1 == 0 || rec.Flag&0x2 != 0) {
			res.regions[readGenicRegion(opt.genes, rec)]++
		}
	}

	return res, nil
}

// javaFormat2 renders v as Java's String.format("%.2f", v) does, which is not
// Go's: Java rounds half-up from the shortest decimal representation of the
// double, so 1.005 and 2.675 become 1.01 and 2.68 (Go gives 1.00 and 2.67), a
// negative that rounds to zero keeps its sign ("-0.00"), and NaN and the
// infinities print as words.
func javaFormat2(v float64) string {
	switch {
	case math.IsNaN(v):
		return "NaN"
	case math.IsInf(v, 1):
		return "Infinity"
	case math.IsInf(v, -1):
		return "-Infinity"
	}
	r, _ := new(big.Rat).SetString(strconv.FormatFloat(math.Abs(v), 'g', -1, 64))
	r.Mul(r, big.NewRat(100, 1))
	r.Add(r, big.NewRat(1, 2))
	cents := new(big.Int).Quo(r.Num(), r.Denom()) // floor, r >= 0
	str := cents.String()
	for len(str) < 3 {
		str = "0" + str
	}
	out := str[:len(str)-2] + "." + str[len(str)-2:]
	if math.Signbit(v) {
		out = "-" + out
	}
	return out
}

// javaStringHash is Java's String.hashCode for an ASCII string.
func javaStringHash(s string) int32 {
	var h int32
	for i := 0; i < len(s); i++ {
		h = 31*h + int32(s[i])
	}
	return h
}

// javaHashOrder returns keys -- given in insertion order -- in the order a
// java.util.HashMap built from them iterates: by bucket, then by insertion
// within a bucket. ngsutilsj prints its tag sections and lists its references
// from HashMaps, so matching its output means matching this.
func javaHashOrder(keys []string) []string {
	capacity := 16
	for float64(len(keys)) > 0.75*float64(capacity) {
		capacity *= 2
	}
	type kb struct {
		key    string
		bucket int
		idx    int
	}
	all := make([]kb, len(keys))
	for i, k := range keys {
		h := uint32(javaStringHash(k))
		h ^= h >> 16
		all[i] = kb{k, int(h & uint32(capacity-1)), i}
	}
	sort.SliceStable(all, func(i, j int) bool { return all[i].bucket < all[j].bucket })
	out := make([]string, len(all))
	for i, e := range all {
		out[i] = e.key
	}
	return out
}

// naturalLess ports compgen-common's natural sort: names split into runs of
// digits and non-digits, digit runs compared as numbers and before any text,
// text runs compared as strings, and a name that runs out first sorting first.
func naturalLess(a, b string) bool {
	ta, tb := naturalTokens(a), naturalTokens(b)
	for i := 0; i < len(ta) || i < len(tb); i++ {
		switch {
		case i >= len(ta):
			return true
		case i >= len(tb):
			return false
		}
		x, y := ta[i], tb[i]
		switch {
		case x.digits && !y.digits:
			return true
		case !x.digits && y.digits:
			return false
		case x.digits:
			if x.num != y.num {
				return x.num < y.num
			}
		default:
			if x.str != y.str {
				return x.str < y.str
			}
		}
	}
	return false
}

type naturalToken struct {
	digits bool
	num    int64
	str    string
}

func naturalTokens(s string) []naturalToken {
	var out []naturalToken
	start := 0
	for i := 1; i <= len(s); i++ {
		if i == len(s) || isDigit(s[i]) != isDigit(s[start]) {
			run := s[start:i]
			if isDigit(s[start]) {
				n, _ := strconv.ParseInt(run, 10, 64)
				out = append(out, naturalToken{digits: true, num: n})
			} else {
				out = append(out, naturalToken{str: run})
			}
			start = i
		}
	}
	return out
}

func isDigit(c byte) bool { return c >= '0' && c <= '9' }

func writeStats(w io.Writer, res *statsResult, opt statsOpts) {
	fmt.Fprintf(w, "Total-reads:\t%d\n", res.total)
	fmt.Fprintf(w, "Mapped-reads:\t%d\n", res.mapped)
	fmt.Fprintf(w, "Unmapped-reads:\t%d\n", res.unmapped)
	fmt.Fprintf(w, "Multiple-mapped-reads:\t%d\n", res.multiple)
	fmt.Fprintf(w, "Uniquely-mapped-reads:\t%d\n", res.mapped-res.multiple)
	fmt.Fprintf(w, "Q30-pct:\t%s%%\n", javaFormat2(100.0*float64(res.q30Bases)/float64(res.totalBases)))
	fmt.Fprintf(w, "Total-bases:\t%d\n", res.totalBases)
	fmt.Fprintf(w, "Ref-length:\t%d\n", res.refLength)
	if !res.hasGaps {
		// With gaps this is RNA-seq, where depth over the genome means nothing.
		fmt.Fprintf(w, "Total-depth:\t%sX\n", javaFormat2(float64(res.totalBases)/float64(res.refLength)))
	}
	if res.paired && res.calcInsertActive {
		fmt.Fprintln(w)
		fmt.Fprintf(w, "Median insert size:\t%d\n", res.inserts.median())
	}

	fmt.Fprintln(w)
	fmt.Fprintln(w, "[Flags]")
	for _, fb := range flagBits {
		if c := res.flagCounts[fb.bit]; c > 0 {
			fmt.Fprintf(w, "%s (0x%x):\t%d\n", fb.name, fb.bit, c)
		}
	}

	for _, name := range javaHashOrder(res.strOrder) {
		fmt.Fprintf(w, "\n[%s]\n", name)
		res.strTags[name].write(w)
	}
	for _, name := range javaHashOrder(res.intOrder) {
		fmt.Fprintf(w, "\n[%s]\n", name)
		res.intTags[name].write(w)
	}

	fmt.Fprintln(w)
	fmt.Fprintln(w, "[References]")
	refs := javaHashOrder(res.refOrder)
	sort.SliceStable(refs, func(i, j int) bool { return naturalLess(refs[i], refs[j]) })
	for _, ref := range refs {
		if c := res.refCounts[ref]; opt.showUnmappedRef || c > 0 {
			fmt.Fprintf(w, "%s\t%d\n", ref, c)
		}
	}
	fmt.Fprintln(w)

	if res.regions != nil {
		writeGeneRegions(w, res.regions)
	}
}

// writeGeneRegions ports the --gtf section of ngsutilsj bam-stats, quirks
// included -- which junction counts each orientation sums, how a zero
// denominator prints -- because pipelines parse it.
func writeGeneRegions(w io.Writer, counts map[gtf.GenicRegion]int64) {
	regions := gtf.GenicRegions()
	maxlen := 0
	for _, r := range regions {
		maxlen = max(maxlen, len(r.Description()))
	}
	fmt.Fprintln(w, "[Gene regions]")
	for _, r := range regions {
		fmt.Fprintf(w, "%-*s\t%d\n", maxlen, r.Description(), counts[r])
	}
	fmt.Fprintln(w)

	sum := func(rs ...gtf.GenicRegion) int64 {
		var n int64
		for _, r := range rs {
			n += counts[r]
		}
		return n
	}
	var sense, antisense int64
	for _, r := range regions {
		if r.IsGene {
			if r.IsSense {
				sense += counts[r]
			} else {
				antisense += counts[r]
			}
		}
	}

	// A 10-fold excess of one strand calls the library stranded. The division
	// is floating point as in Java: no antisense reads gives +Inf (FR), none of
	// either gives NaN, which matches neither test and so reads as unstranded.
	var intronic, exonic, junction int64
	strandedness := math.Log10(float64(sense) / float64(antisense))
	switch {
	case strandedness > 1:
		intronic = sum(gtf.NCIntron, gtf.CodingIntron, gtf.UTR3Intron, gtf.UTR5Intron)
		exonic = sense - intronic
		junction = sum(gtf.Junction, gtf.JunctionAnti)
		fmt.Fprintln(w, "Orientation\tFR")
	case strandedness < -1:
		intronic = sum(gtf.NCIntronAnti, gtf.CodingIntronAnti, gtf.UTR3IntronAnti, gtf.UTR5IntronAnti)
		exonic = antisense - intronic
		junction = sum(gtf.JunctionAnti, gtf.NCJunctionAnti)
		fmt.Fprintln(w, "Orientation\tRF")
	default:
		intronic = sum(gtf.NCIntron, gtf.CodingIntron, gtf.UTR3Intron, gtf.UTR5Intron,
			gtf.NCIntronAnti, gtf.CodingIntronAnti, gtf.UTR3IntronAnti, gtf.UTR5IntronAnti)
		exonic = sense + antisense - intronic
		junction = sum(gtf.Junction, gtf.JunctionAnti)
		fmt.Fprintln(w, "Orientation\tunstranded")
	}

	fmt.Fprintf(w, "Sense      \t%d\n", sense)
	fmt.Fprintf(w, "Anti-sense \t%d\n", antisense)
	fmt.Fprintln(w)
	fmt.Fprintf(w, "Exonic     \t%d\n", exonic)
	fmt.Fprintf(w, "Intronic   \t%d\n", intronic)
	fmt.Fprintln(w)
	fmt.Fprintf(w, "Junction   \t%d\n", junction)
	fmt.Fprintln(w)

	if antisense > 0 {
		ratio := -float64(antisense) / float64(sense)
		if sense > antisense {
			ratio = float64(sense) / float64(antisense)
		}
		fmt.Fprintf(w, "Sense/anti-sense ratio     \t%s\n", javaFormat2(ratio))
	} else {
		fmt.Fprintln(w, "Sense/anti-sense ratio     \t0")
	}
	if intronic > 0 {
		fmt.Fprintf(w, "Exonic/intronic ratio      \t%s\n", javaFormat2(float64(exonic)/float64(intronic)))
	} else {
		fmt.Fprintln(w, "Exonic/intronic ratio      \t0")
	}
	if ig := counts[gtf.Intergenic]; ig > 0 {
		fmt.Fprintf(w, "Exonic/genomic ratio       \t%s\n", javaFormat2(float64(exonic)/float64(ig)))
	} else {
		fmt.Fprintln(w, "Exonic/genomic ratio       \t0")
	}
	if junction > 0 {
		fmt.Fprintf(w, "Non-junction/junction ratio\t%s\n", javaFormat2(float64(exonic-junction)/float64(junction)))
	} else {
		fmt.Fprintln(w, "Non-junction/junction ratio\t0")
	}
}
