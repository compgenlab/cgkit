package samcmd

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/htsio"
	"github.com/compgenlab/cghts/vcf"
	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/compgenlab/cgkit/internal/cmdio"
	"github.com/spf13/cobra"
)

// defaultVcfCountFlagExclude drops unmapped (0x4), secondary (0x100), QC-fail
// (0x200), duplicate (0x400) and supplementary (0x800) reads -- the same 3844
// ngsutilsj's vcf-count defaulted to. A secondary or supplementary record is a
// second copy of a read already counted at its primary location, and counting
// duplicates inflates whichever allele the duplicated fragment carried.
const defaultVcfCountFlagExclude = 0x4 | 0x100 | 0x200 | 0x400 | 0x800

type vcfCountOpts struct {
	output       string
	flagRequired int
	flagExclude  int
	minMapQ      int
	minBaseQual  int
	properPairs  bool
	bothMates    bool
	cramRef      string
	tags         tagFilterFlags

	region      string
	passing     bool
	snvsOnly    bool
	indelsOnly  bool
	het         bool
	sample      string
	skipMissing bool
	batchLen    int

	outID    bool
	outAF    bool
	outTotal bool
	outOther bool
	outGT    bool
	verbose  bool
}

var vcfCount vcfCountOpts

var samVcfCountCmd = &cobra.Command{
	GroupID:     "samcmd",
	Annotations: map[string]string{"since": "v0.6.11"},
	Use:         "sam-vcfcount <input.bam> <input.vcf>",
	Short:       "Count reference and alternate reads for each VCF variant in a SAM/BAM/CRAM",
	Long: `For each variant in a VCF, count the reads in a SAM/BAM/CRAM file that carry
the reference allele and the alternate allele. One row is written per ALT
allele, so a multiallelic record yields one row per alternate, each with the
same ref_count.

Reads are classified whole rather than base-by-base: the bases a read aligns
across the record's REF span (plus any insertion inside it) are compared to REF
and to each ALT, so SNVs, MNVs, insertions, deletions and complex substitutions
are all matched exactly. A read that spans the variant and matches neither is
"other" (--other reports it); a read that does not span it is not counted.
Indels are matched as the aligner placed them: no realignment is done, so an
indel the aligner put elsewhere in a repeat is not recognised.

A read is counted only if it passes the read filters (--flag-exclude,
--min-mapq, tag filters, ...) and, for SNVs/MNVs, every compared base has
quality >= --min-basequal. Base quality is not applied to indels, matching
samtools mpileup. By default a fragment is counted once even when both mates
overlap the variant; mates that disagree are not counted. --both-mates counts
each read separately.

The alignment file must be coordinate-sorted and indexed (.bai/.crai). The VCF
may be streamed (or "-" for stdin); its contig names are matched to the
alignment file's whether they are spelled chr1, 1 or NC_000001.11. A variant
with no reads gets counts of 0.

Not counted: symbolic alleles (<DEL>, <NON_REF>, ...), breakends, the spanning
deletion "*", and gVCF reference blocks.

Output columns: chrom pos [id] ref alt ref_count alt_count [alt_freq]
[total_count] [other_count] [gt]. alt_freq is alt_count over every read that
spanned the variant (total_count), so at a multiallelic site it is not
alt/(alt+ref).`,
	Args: cobra.ExactArgs(2),
	RunE: func(cmd *cobra.Command, args []string) error {
		return runSamVcfCount(cmd, args[0], args[1], &vcfCount)
	},
}

func init() {
	f := samVcfCountCmd.Flags()
	o := &vcfCount
	f.StringVarP(&o.output, "output", "o", "", "Output file (default stdout)")

	// Read filters.
	f.IntVar(&o.flagRequired, "flag-required", 0, "Count reads with all of these flag bits set")
	f.IntVar(&o.flagExclude, "flag-exclude", defaultVcfCountFlagExclude, "Exclude reads with any of these flag bits set (default excludes unmapped, secondary, QC-fail, duplicate, supplementary)")
	f.IntVar(&o.minMapQ, "min-mapq", 0, "Count reads with mapping quality at or above this value")
	f.IntVar(&o.minBaseQual, "min-basequal", 13, "Minimum base quality for the bases compared at an SNV/MNV")
	f.BoolVar(&o.properPairs, "proper-pairs", false, "Only count properly paired reads (flag 0x2 set, mate mapped)")
	f.BoolVar(&o.bothMates, "both-mates", false, "Count both mates of an overlapping pair rather than the fragment once")
	f.StringVar(&o.cramRef, "cram-ref", "", "Reference FASTA for CRAM files (path, http(s):// URL, or s3://)")
	o.tags.register(samVcfCountCmd)

	// Variant filters.
	f.StringVar(&o.region, "region", "", "Only count variants in this region (chrom, or chrom:start-end, 1-based)")
	f.BoolVar(&o.passing, "passing", false, "Only count variants with FILTER PASS or '.'")
	f.BoolVar(&o.snvsOnly, "snvs", false, "Only count SNVs")
	f.BoolVar(&o.indelsOnly, "indels", false, "Only count indels")
	f.BoolVar(&o.het, "het", false, "Only count variants heterozygous in --sample (requires GT)")
	f.StringVar(&o.sample, "sample", "", "VCF sample for --het and --gt (default: the first sample)")
	f.BoolVar(&o.skipMissing, "skip-missing", false, "Skip variants on contigs absent from the alignment file (default: error)")
	f.IntVar(&o.batchLen, "batch-len", 10000, "Query the alignment file once per window of up to this many bases")

	// Extra output columns.
	f.BoolVar(&o.outID, "id", false, "Add the VCF ID column")
	f.BoolVar(&o.outAF, "af", false, "Add an alt_freq column (alt_count / total_count)")
	f.BoolVar(&o.outTotal, "total", false, "Add a total_count column (reads spanning the variant)")
	f.BoolVar(&o.outOther, "other", false, "Add an other_count column (spanning reads matching neither REF nor this ALT)")
	f.BoolVar(&o.outGT, "gt", false, "Add the --sample genotype from the VCF")
	f.BoolVarP(&o.verbose, "verbose", "v", false, "Report a summary of what was counted and skipped on stderr")
	samVcfCountCmd.MarkFlagsMutuallyExclusive("snvs", "indels")
}

// vcVariant is one VCF record prepared for counting.
type vcVariant struct {
	rec     *vcf.VcfRecord
	bamRef  string   // contig name as the alignment file spells it
	start   int      // 0-based start of REF
	end     int      // 0-based exclusive end of REF
	ref     string   // upper-cased REF
	alts    []string // upper-cased ALTs; "" for an allele that cannot be counted
	isIndel bool     // any countable ALT differs in length from REF
	gt      string

	// Results, filled by countVariants.
	refCount   int
	altCounts  []int
	total      int            // reads that spanned the variant, whatever they carried
	matesFound map[string]int // read name -> allele index (-1 ref, -2 other, -3 conflict)
}

// isCountableAlt reports whether alt is a literal allele whose reads can be
// matched: not symbolic, not a breakend, not the spanning deletion or missing.
func isCountableAlt(alt string) bool {
	if alt == "" || alt == "." || alt == "*" {
		return false
	}
	for i := 0; i < len(alt); i++ {
		switch alt[i] {
		case 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n':
		default:
			return false
		}
	}
	return true
}

type vcfCountSummary struct {
	records, counted, filtered, uncountable, missingContig int
	contigsMissing                                         map[string]bool
}

func runSamVcfCount(cmd *cobra.Command, bamPath, vcfPath string, o *vcfCountOpts) error {
	if o.batchLen < 0 {
		return fmt.Errorf("--batch-len must be >= 0")
	}
	tagFilters, err := o.tags.parse()
	if err != nil {
		return err
	}
	opts := htsio.NewSamReaderOpts()
	flagReq, flagEx := o.flagRequired, o.flagExclude
	if o.properPairs {
		flagReq |= 0x2
		flagEx |= 0x8
	}
	if flagReq != 0 {
		opts.FlagRequired(flagReq)
	}
	if flagEx != 0 {
		opts.FlagFilter(flagEx)
	}
	if o.minMapQ > 0 {
		opts.MinMapQ(o.minMapQ)
	}
	for _, tf := range tagFilters {
		opts.AddTagFilter(tf)
	}
	if o.cramRef != "" {
		opts.RefPath(o.cramRef)
	}

	var regionRef string
	regionStart, regionEnd := 0, -1
	if o.region != "" {
		regionRef, regionStart, regionEnd, err = htsio.ParseRegion(o.region)
		if err != nil {
			return err
		}
	}

	reader, err := htsio.OpenSamReader(cmd.Context(), bamPath, opts)
	if err != nil {
		return fmt.Errorf("opening %s: %w", bamPath, err)
	}
	defer reader.Close()
	samHeader, err := reader.Header()
	if err != nil {
		return fmt.Errorf("reading header of %s: %w", bamPath, err)
	}
	var bamRefs []string
	for _, r := range samHeader.References() {
		bamRefs = append(bamRefs, r.Name)
	}
	contigs := vcf.NewContigConverter(bamRefs)

	var vr *vcf.VcfReader
	if vcfPath == "-" {
		vr, err = vcf.NewVcfReader(cmd.InOrStdin())
	} else {
		vr, err = vcf.OpenVcfFile(cmd.Context(), vcfPath)
	}
	if err != nil {
		return fmt.Errorf("opening %s: %w", vcfPath, err)
	}
	defer vr.Close()
	vh, err := vr.Header()
	if err != nil {
		return fmt.Errorf("reading header of %s: %w", vcfPath, err)
	}
	sampleIdx := -1
	if o.het || o.outGT {
		if len(vh.Samples()) == 0 {
			return fmt.Errorf("--het and --gt need a VCF with samples; %s has none", vcfPath)
		}
		sampleIdx = 0
		if o.sample != "" {
			sampleIdx = vh.SampleIndex(o.sample)
			if sampleIdx < 0 {
				return fmt.Errorf("--sample %q is not in %s (samples: %s)", o.sample, vcfPath, strings.Join(vh.Samples(), ", "))
			}
		}
	} else if o.sample != "" {
		return fmt.Errorf("--sample only applies with --het or --gt")
	}

	out, err := cmdio.CreateOutput(cmd, "-o/--output", o.output)
	if err != nil {
		return err
	}
	w := bufio.NewWriter(out.W)
	fail := func(err error) error {
		out.Discard()
		return err
	}

	fmt.Fprintln(w, "## program: "+buildinfo.String())
	fmt.Fprintln(w, "## cmd: "+buildinfo.CommandLine())
	fmt.Fprintln(w, "## vcf-input: "+vcfPath)
	fmt.Fprintln(w, "## bam-input: "+bamPath)
	cols := []string{"chrom", "pos"}
	if o.outID {
		cols = append(cols, "id")
	}
	cols = append(cols, "ref", "alt", "ref_count", "alt_count")
	if o.outAF {
		cols = append(cols, "alt_freq")
	}
	if o.outTotal {
		cols = append(cols, "total_count")
	}
	if o.outOther {
		cols = append(cols, "other_count")
	}
	if o.outGT {
		cols = append(cols, "gt")
	}
	fmt.Fprintln(w, strings.Join(cols, "\t"))

	sum := vcfCountSummary{contigsMissing: map[string]bool{}}
	var batch []*vcVariant
	flush := func() error {
		if len(batch) == 0 {
			return nil
		}
		if err := countBatch(reader, batch, o.minBaseQual, !o.bothMates); err != nil {
			return err
		}
		for _, v := range batch {
			writeVcfCountRows(w, v, o)
		}
		batch = batch[:0]
		return nil
	}

	for {
		rec, err := vr.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fail(err)
		}
		sum.records++
		if o.region != "" {
			if !sameContig(rec.Chrom, regionRef) || rec.Pos <= regionStart || (regionEnd >= 0 && rec.Pos > regionEnd) {
				continue
			}
		}
		if o.passing && rec.IsFiltered() {
			sum.filtered++
			continue
		}
		v := newVcVariant(rec)
		if v == nil {
			sum.uncountable++
			continue
		}
		if (o.snvsOnly && v.isIndel) || (o.indelsOnly && !v.isIndel) {
			sum.filtered++
			continue
		}
		if sampleIdx >= 0 {
			gt, err := sampleGT(rec, sampleIdx)
			if err != nil {
				return fail(err)
			}
			if o.het {
				if gt == "" {
					return fail(fmt.Errorf("--het: no GT for the sample at %s:%d", rec.Chrom, rec.Pos))
				}
				if !isHet(gt) {
					sum.filtered++
					continue
				}
			}
			v.gt = gt
		}
		bamRef, ok := contigs.Resolve(rec.Chrom)
		if !ok {
			if !o.skipMissing {
				return fail(fmt.Errorf("contig %q (%s:%d) is not in %s; use --skip-missing to skip such variants", rec.Chrom, rec.Chrom, rec.Pos, bamPath))
			}
			sum.missingContig++
			sum.contigsMissing[rec.Chrom] = true
			continue
		}
		v.bamRef = bamRef

		// Start a new window on a contig change, a step backwards, or once the
		// window would exceed --batch-len; an unsorted VCF just costs queries.
		if len(batch) > 0 {
			first := batch[0]
			if first.bamRef != v.bamRef || v.start < first.start || v.start-first.start > o.batchLen {
				if err := flush(); err != nil {
					return fail(err)
				}
			}
		}
		batch = append(batch, v)
		sum.counted++
	}
	if err := flush(); err != nil {
		return fail(err)
	}
	if err := w.Flush(); err != nil {
		return fail(err)
	}
	if err := out.Close(); err != nil {
		return err
	}

	if o.verbose {
		e := cmd.ErrOrStderr()
		fmt.Fprintf(e, "sam-vcfcount: %d records read, %d counted, %d filtered out, %d with no countable allele\n",
			sum.records, sum.counted, sum.filtered, sum.uncountable)
		if sum.missingContig > 0 {
			var names []string
			for c := range sum.contigsMissing {
				names = append(names, c)
			}
			fmt.Fprintf(e, "sam-vcfcount: %d variants skipped on contigs absent from %s: %s\n",
				sum.missingContig, bamPath, strings.Join(names, ", "))
		}
	}
	return nil
}

func sameContig(a, b string) bool {
	if a == b {
		return true
	}
	ka, oka := vcf.CanonicalContig(a)
	kb, okb := vcf.CanonicalContig(b)
	return oka && okb && ka == kb
}

// newVcVariant prepares rec for counting, or returns nil when no ALT allele can
// be matched against reads.
func newVcVariant(rec *vcf.VcfRecord) *vcVariant {
	if rec.IsRefBlock() || !isCountableAlt(rec.Ref) {
		return nil
	}
	v := &vcVariant{
		rec:   rec,
		start: rec.Pos - 1,
		end:   rec.Pos - 1 + len(rec.Ref),
		ref:   strings.ToUpper(rec.Ref),
	}
	any := false
	for _, a := range rec.Alt() {
		if !isCountableAlt(a) {
			v.alts = append(v.alts, "")
			continue
		}
		any = true
		a = strings.ToUpper(a)
		if len(a) != len(v.ref) {
			v.isIndel = true
		}
		v.alts = append(v.alts, a)
	}
	if !any {
		return nil
	}
	v.altCounts = make([]int, len(v.alts))
	return v
}

func sampleGT(rec *vcf.VcfRecord, idx int) (string, error) {
	s, err := rec.Sample(idx)
	if err != nil {
		return "", err
	}
	gt, ok := s.Get("GT")
	if !ok || gt.IsMissing() {
		return "", nil
	}
	return gt.String(), nil
}

// isHet reports whether gt names at least two different called alleles.
func isHet(gt string) bool {
	alleles := strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
	first := ""
	for _, a := range alleles {
		if a == "." {
			continue
		}
		if first == "" {
			first = a
		} else if a != first {
			return true
		}
	}
	return false
}

func writeVcfCountRows(w io.Writer, v *vcVariant, o *vcfCountOpts) {
	for i, alt := range v.alts {
		if alt == "" {
			continue
		}
		other := v.total - v.refCount - v.altCounts[i]
		row := []string{v.rec.Chrom, strconv.Itoa(v.rec.Pos)}
		if o.outID {
			id := v.rec.ID()
			if id == "" {
				id = "."
			}
			row = append(row, id)
		}
		row = append(row, v.rec.Ref, v.rec.Alt()[i], strconv.Itoa(v.refCount), strconv.Itoa(v.altCounts[i]))
		if o.outAF {
			if v.total > 0 {
				row = append(row, strconv.FormatFloat(float64(v.altCounts[i])/float64(v.total), 'g', 6, 64))
			} else {
				row = append(row, ".")
			}
		}
		if o.outTotal {
			row = append(row, strconv.Itoa(v.total))
		}
		if o.outOther {
			row = append(row, strconv.Itoa(other))
		}
		if o.outGT {
			gt := v.gt
			if gt == "" {
				gt = "."
			}
			row = append(row, gt)
		}
		fmt.Fprintln(w, strings.Join(row, "\t"))
	}
}

// countBatch queries the reads over a window of variants on one contig and
// tallies each variant. Variants must share bamRef.
func countBatch(reader htsio.SamReader, batch []*vcVariant, minBQ int, oncePerFragment bool) error {
	qStart, qEnd := batch[0].start, batch[0].end
	for _, v := range batch[1:] {
		qStart = min(qStart, v.start)
		qEnd = max(qEnd, v.end)
	}
	// One base either side: an insertion after the last REF base is only seen
	// by a read reaching past it.
	if qStart > 0 {
		qStart--
	}
	qEnd++
	recs, err := reader.Query(batch[0].bamRef, qStart, qEnd)
	if err != nil {
		return fmt.Errorf("querying %s:%d-%d (is the alignment file indexed?): %w", batch[0].bamRef, qStart+1, qEnd, err)
	}
	var reads []*alignedRead
	for rec, err := range recs {
		if err != nil {
			return err
		}
		if ar := newAlignedRead(rec); ar != nil {
			reads = append(reads, ar)
		}
	}
	countVariants(batch, reads, minBQ, oncePerFragment)
	return nil
}

// countVariants tallies reads against each variant.
func countVariants(batch []*vcVariant, reads []*alignedRead, minBQ int, oncePerFragment bool) {
	for _, v := range batch {
		if oncePerFragment {
			v.matesFound = map[string]int{}
		}
		for _, r := range reads {
			if r.start > v.start {
				break // reads arrive sorted by start; none later can span v
			}
			if r.end <= v.start {
				continue
			}
			allele, ok := r.alleleAt(v, minBQ)
			if !ok {
				continue
			}
			idx := -2 // other
			if allele == v.ref {
				idx = -1
			} else {
				for i, a := range v.alts {
					if a != "" && a == allele {
						idx = i
						break
					}
				}
			}
			if oncePerFragment && r.paired {
				if prev, seen := v.matesFound[r.name]; seen {
					if prev != idx && prev != -3 {
						v.matesFound[r.name] = -3 // conflict: neither mate counts
					}
					continue
				}
				v.matesFound[r.name] = idx
				continue
			}
			v.tally(idx)
		}
		for _, idx := range v.matesFound {
			if idx != -3 {
				v.tally(idx)
			}
		}
		v.matesFound = nil
	}
}

func (v *vcVariant) tally(idx int) {
	v.total++
	switch {
	case idx == -1:
		v.refCount++
	case idx >= 0:
		v.altCounts[idx]++
	}
}

// alignedRead is a read reduced to what allele matching needs.
type alignedRead struct {
	name   string
	paired bool
	start  int // 0-based first aligned reference base
	end    int // 0-based exclusive
	seq    string
	qual   string // "" when absent
	ops    []htsio.CigarOp
}

func newAlignedRead(rec *htsio.SamRecord) *alignedRead {
	if rec.IsUnmapped() || rec.Seq == "" || rec.Seq == "*" {
		return nil
	}
	ops, err := htsio.ParseCigar(rec.Cigar)
	if err != nil || len(ops) == 0 {
		return nil
	}
	q := rec.Qual
	if q == "*" || len(q) != len(rec.Seq) {
		q = ""
	}
	return &alignedRead{
		name:   rec.ReadName,
		paired: rec.Flag&0x1 != 0,
		start:  rec.AlignStart(),
		end:    rec.AlignEnd(),
		seq:    strings.ToUpper(rec.Seq),
		qual:   q,
		ops:    ops,
	}
}

// alleleAt returns the sequence this read carries over v's REF span: the bases
// aligned to [start,end) plus any insertion that falls inside it or right
// after its last base. A read deleting the whole span carries "", which
// matches no allele and counts as other. ok is false when the read does not
// span the variant (it starts or ends inside it, or is spliced across it) or,
// for an SNV/MNV, a compared base is below minBQ.
//
// An indel's read must also reach the base after REF, or a read ending on an
// insertion's anchor would read as REF. Both flanks -- the first REF base,
// which VCF uses as the anchor, and the base after REF -- must be aligned
// rather than deleted, or a read with a longer deletion than the variant's
// would read as ALT; such a read is other.
func (r *alignedRead) alleleAt(v *vcVariant, minBQ int) (string, bool) {
	right := v.end - 1
	if v.isIndel {
		right = v.end
	}
	if r.start > v.start || r.end <= right {
		return "", false
	}

	var b strings.Builder
	refPos, qPos := r.start, 0
	leftAligned, rightAligned := !v.isIndel, !v.isIndel
	checkBQ := !v.isIndel && r.qual != "" && minBQ > 0
	for _, op := range r.ops {
		if refPos > right {
			break
		}
		switch op.Op {
		case 'M', '=', 'X':
			for i := 0; i < op.Len; i++ {
				p := refPos + i
				if p >= v.start && p < v.end {
					if checkBQ && int(r.qual[qPos+i])-33 < minBQ {
						return "", false
					}
					b.WriteByte(r.seq[qPos+i])
				}
				if p == v.start {
					leftAligned = true
				}
				if p == right {
					rightAligned = true
				}
			}
			refPos += op.Len
			qPos += op.Len
		case 'I':
			// An insertion sits between refPos-1 and refPos.
			if refPos-1 >= v.start && refPos-1 < v.end {
				b.WriteString(r.seq[qPos : qPos+op.Len])
			}
			qPos += op.Len
		case 'S':
			qPos += op.Len
		case 'D':
			// Deleted bases inside the span contribute nothing.
			refPos += op.Len
		case 'N':
			if refPos <= right && refPos+op.Len > v.start {
				return "", false // spliced across the variant
			}
			refPos += op.Len
		}
	}
	if !leftAligned || !rightAligned {
		// Spans the variant but deletes a flank: a different event, so it
		// carries neither allele.
		return "", true
	}
	return b.String(), true
}
