package ontcmd

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/align"
	"github.com/compgenlab/cghts/htsio"
	_ "github.com/compgenlab/cghts/htsio/bam"
	_ "github.com/compgenlab/cghts/htsio/cram"
	_ "github.com/compgenlab/cghts/htsio/sam"
	"github.com/compgenlab/cghts/seqio"
	"github.com/compgenlab/cghts/support/sequtils"
	"github.com/spf13/cobra"
)

// polyaGeom is the direction-dependent geometry of a poly(A) search, in SEQ
// (reference-forward) coordinates. SEQ index always increases with reference
// coordinate, so strand enters the search only through these three fields.
type polyaGeom struct {
	tailBase byte // 'A' when the mRNA is on the + strand, 'T' when on the -
	terminus int  // SEQ index of the mRNA 3'-most base
	step     int  // SEQ index delta when walking 5'-ward along the mRNA
}

// geomFor returns the search geometry for a read of seqLen bases whose mRNA is
// on the plus (mRNAPlus) or minus strand.
//
// A minus-strand mRNA means BAM stored revcomp(read), so the tail sits at the
// start of SEQ as T bases and mRNA 5'->3' runs toward decreasing SEQ index.
func geomFor(seqLen int, mRNAPlus bool) polyaGeom {
	if mRNAPlus {
		return polyaGeom{tailBase: 'A', terminus: seqLen - 1, step: -1}
	}
	return polyaGeom{tailBase: 'T', terminus: 0, step: 1}
}

// polyaParams are the windowed A-fraction detection criteria.
type polyaParams struct {
	window   int
	minAFrac float64
	minLen   int
	maxJunk  int
}

// polyaCall is a detector result, in SEQ coordinates. start is the first tail
// base in mRNA 5'->3' order, so start > end when the mRNA is on the minus strand.
type polyaCall struct {
	found  bool
	start  int
	end    int
	length int
}

func absInt(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// isTailBase reports whether c is the tail base, folding case.
//
// This is a plain byte compare rather than sequutils.DNAMatches on purpose:
// DNAMatches is IUPAC-aware, so DNAMatches('N', 'A') is true, and N bases would
// inflate the A-fraction and pull the trace through ambiguous sequence.
func isTailBase(c, tailBase byte) bool {
	if c >= 'a' && c <= 'z' {
		c -= 'a' - 'A'
	}
	return c == tailBase
}

// windowAFrac returns the fraction of tail bases in the window of `window` bases
// running from SEQ index i back toward end.
//
// The window looks backward, toward the anchor, never ahead: a lookahead window
// starts failing a full window's width before the real tail boundary.
//
// The denominator is always `window`, even near end where fewer bases are
// available. Positions beyond the anchor -- off the read end, or 3'-ward junk
// the anchor scan already rejected -- are counted as tail rather than as misses,
// because they are not evidence against the tail. Without this a single miscall
// near the 3' tip fails a 2- or 3-base partial window (one T in four bases is
// only 0.75) and costs the anchor, so the retry re-anchors inward and the
// reported length loses the tip.
func windowAFrac(seq string, g polyaGeom, window, i, end int) float64 {
	if window <= 0 {
		return 0
	}
	avail := absInt(i-end) + 1
	if avail > window {
		avail = window
	}
	hits := 0
	j := i
	for k := 0; k < avail; k++ {
		if isTailBase(seq[j], g.tailBase) {
			hits++
		}
		j -= g.step // back toward end
	}
	hits += window - avail // positions beyond the anchor: not misses
	return float64(hits) / float64(window)
}

// extendTail walks 5'-ward from anchor while the trailing window keeps its
// A-fraction, and returns the SEQ index of the last base that is itself a tail
// base and whose window passed.
//
// Walking continues through the soft-clip boundary into aligned bases, so
// genome-encoded A's that the aligner absorbed into the alignment are counted as
// part of the tail. Anchoring the result on a real tail base rather than
// wherever the window happened to fail is what keeps a large --window from
// inflating the reported length: the walk may overrun into genomic sequence, but
// the result rewinds to the last tail base.
func extendTail(seq string, g polyaGeom, p polyaParams, anchor int) int {
	start := anchor // the caller guarantees anchor is itself a tail base
	for i := anchor + g.step; i >= 0 && i < len(seq); i += g.step {
		if windowAFrac(seq, g, p.window, i, anchor) < p.minAFrac {
			break
		}
		if isTailBase(seq[i], g.tailBase) {
			start = i
		}
	}
	return start
}

// findPolyA locates the poly(A) tail in seq (reference-forward SEQ, as stored in
// the BAM). scanFrom is where the junk scan begins: the read terminus, or the
// base just inside a matched 3' adapter.
//
// Candidate anchors are tried outermost-first, and an anchor whose extension
// falls short of p.minLen is retried from the next tail base inward until the
// p.maxJunk budget is spent. The retry is not optional: an untrimmed adapter
// ending in a lone A (say "...[A x 20][GGATCA]") would otherwise anchor on the
// adapter's A, fail the window immediately on the C, and report no call — the
// real 20 nt tail never examined. The retry cannot invent a tail, since reaching
// minLen from a spurious anchor requires minLen A-rich bases inward of it.
func findPolyA(seq string, g polyaGeom, p polyaParams, scanFrom int) polyaCall {
	if len(seq) == 0 || seq == "*" {
		return polyaCall{}
	}
	if scanFrom < 0 || scanFrom >= len(seq) {
		return polyaCall{}
	}
	junk := 0
	for i := scanFrom; i >= 0 && i < len(seq); {
		if !isTailBase(seq[i], g.tailBase) {
			junk++
			if junk > p.maxJunk {
				return polyaCall{}
			}
			i += g.step
			continue
		}
		start := extendTail(seq, g, p, i)
		length := absInt(i-start) + 1
		if length >= p.minLen {
			return polyaCall{found: true, start: start, end: i, length: length}
		}
		// Too short to be a tail: charge the whole candidate run to the junk
		// budget and resume scanning inward from just past it.
		junk += length
		if junk > p.maxJunk {
			return polyaCall{}
		}
		i = start + g.step
	}
	return polyaCall{}
}

// polyaRefPos maps the SEQ index of the first tail base to a 1-based reference
// position, and reports which part of the alignment it fell in.
//
// No strand branching is needed: SEQ index always increases with reference
// coordinate, so the same three cases cover both strands. A tail base in the
// leading clip is genuinely before the alignment, and one in the trailing clip
// is genuinely past it, whichever strand the mRNA is on.
func polyaRefPos(rec *htsio.SamRecord, startIdx, step int) (pos int, src string, ok bool) {
	if rec.Cigar == "*" || rec.Seq == "*" {
		return 0, "", false
	}
	if r0, ok := htsio.CigarQueryToRef(rec.Pos, rec.Cigar, startIdx); ok {
		return r0 + 1, "aligned", true
	}
	lead, trail := htsio.CigarSoftClips(rec.Cigar)
	switch {
	case startIdx < lead:
		// Leading clip: the site is bounded by the first aligned base.
		return rec.Pos - 1, "clip", true
	case startIdx >= len(rec.Seq)-trail:
		// Trailing clip: the site is bounded by the first base past the
		// alignment. When junk sits between the alignment end and the tail those
		// bases are adapter, not genome, so they carry no genomic offset to add.
		return htsio.CigarAlignEnd(rec.Pos, rec.Cigar) + 1, "clip", true
	default:
		// Insertion: anchor on the nearest reference base 5'-ward along the
		// mRNA, i.e. toward the transcript body. Inherently +/-1.
		for j := startIdx + step; j >= 0 && j < len(rec.Seq); j += step {
			if r0, ok := htsio.CigarQueryToRef(rec.Pos, rec.Cigar, j); ok {
				return r0 + 1, "ins", true
			}
		}
		return 0, "", false
	}
}

// terminalHardClip reports whether the CIGAR's outermost operation on the mRNA
// 3' side is a hard clip, meaning the tail was physically removed from SEQ and
// anything found in its place is an artifact.
func terminalHardClip(cigar string, mRNAPlus bool) bool {
	ops, err := htsio.ParseCigar(cigar)
	if err != nil || len(ops) == 0 {
		return false
	}
	if mRNAPlus {
		return ops[len(ops)-1].Op == 'H'
	}
	return ops[0].Op == 'H'
}

// adapterScanStart locates the 3' adapter in a bounded probe window at the read
// terminus and returns the SEQ index just inside it, or g.terminus when no
// adapter is configured or none matched.
func adapterScanStart(seq string, g polyaGeom, aligner align.PairwiseAligner, adapter string, minIdent float64, maxJunk int) int {
	if adapter == "" || aligner == nil || len(seq) == 0 {
		return g.terminus
	}
	// The clip holds [tail][adapter] when the mRNA is on the plus strand and
	// [revcomp(adapter)][tail] when it is on the minus strand.
	probe := adapter
	if !isTailBase('A', g.tailBase) {
		probe = sequtils.ReverseComplement(adapter)
	}
	// Bound the target. ONT reads reach 100 kb and Smith-Waterman is
	// O(|probe| x |target|) per read, so never align against the whole SEQ.
	span := 2*len(adapter) + maxJunk
	lo, hi := 0, len(seq)
	if g.step < 0 {
		if lo = len(seq) - span; lo < 0 {
			lo = 0
		}
	} else {
		if hi = span; hi > len(seq) {
			hi = len(seq)
		}
	}
	target := seq[lo:hi]
	if len(target) == 0 {
		return g.terminus
	}

	// Empty quality strings are safe here: the aligner and Matches() never read
	// qual. Do not call SeqQual.Sub() on these — it slices seq and qual alike and
	// would panic.
	aln := aligner.Align(
		seqio.NewStringSeqQual(probe, "", "adapter").FullSeq(),
		seqio.NewStringSeqQual(target, "", "probe").FullSeq(),
	)
	if aln == nil || float64(aln.Matches())/float64(len(probe)) < minIdent {
		return g.terminus
	}

	// TargetStart/TargetEnd are 0-based half-open within `target`.
	var inner int // the adapter base closest to the tail
	if g.step < 0 {
		inner = lo + aln.TargetStart
	} else {
		inner = lo + aln.TargetEnd - 1
	}
	scan := inner + g.step
	if scan < 0 || scan >= len(seq) {
		return g.terminus
	}
	return scan
}

type polyaStats struct {
	total   int
	called  int
	noCall  int
	skipped int
	clamped int
}

func runPolya(inputFile string) error {
	params := polyaParams{
		window:   polyaWindow,
		minAFrac: polyaMinAFrac,
		minLen:   polyaMinLen,
		maxJunk:  polyaMaxJunk,
	}

	opts := htsio.NewSamReaderOpts()
	if polyaCramRef != "" {
		opts.RefPath(polyaCramRef)
	}
	if polyaThreads > 1 {
		opts.Threads(polyaThreads)
	}
	reader, err := htsio.NewSamReader(inputFile, opts)
	if err != nil {
		return err
	}
	defer reader.Close()

	header, err := reader.Header()
	if err != nil {
		return fmt.Errorf("reading header: %w", err)
	}
	refLens := make(map[string]int)
	if header != nil {
		for _, r := range header.References() {
			refLens[r.Name] = r.Length
		}
	}

	var aligner align.PairwiseAligner
	if polyaAdapter != "" {
		// Built once: the aligner is thread-safe and stateless per Align() call.
		aligner = align.NewLocalAligner(align.OntAlignmentDefaults().ClippingDisable())
	}

	out, closeFn, err := openWriter(polyaOutput, true)
	if err != nil {
		return err
	}
	// One row per read means one write per read, so buffer: ont-umi-lookup gets
	// away without this only because it emits rows for a subset of reads.
	w := bufio.NewWriterSize(out, 64*1024)

	cols := []string{"read_name", "chrom", "polya_pos", "strand"}
	if polyaShowLength {
		cols = append(cols, "polya_len")
	}
	if polyaShowSource {
		cols = append(cols, "polya_source")
	}
	cols = append(cols, polyaTags...)
	fmt.Fprintln(w, strings.Join(cols, "\t"))

	var stats polyaStats
	row := make([]string, 0, len(cols))

	emit := func(rec *htsio.SamRecord, chrom, posStr, strand, lenStr, srcStr string) {
		row = row[:0]
		row = append(row, rec.ReadName, chrom, posStr, strand)
		if polyaShowLength {
			row = append(row, lenStr)
		}
		if polyaShowSource {
			row = append(row, srcStr)
		}
		for _, tg := range polyaTags {
			if t, ok := rec.Tags[tg]; ok {
				row = append(row, t.Value)
			} else {
				row = append(row, "NA")
			}
		}
		fmt.Fprintln(w, strings.Join(row, "\t"))
	}

	for rec, err := range reader.Records() {
		if err != nil {
			closeFn()
			return err
		}
		// Every read has exactly one primary alignment, so skipping these keeps
		// one row per read rather than duplicating read names. Supplementary
		// alignments are usually hard-clipped anyway, so the tail is not there.
		if rec.IsSecondary() || rec.IsSupplementary() {
			stats.skipped++
			continue
		}
		stats.total++

		chrom := rec.RefName
		if chrom == "" || chrom == "*" {
			chrom = "NA"
		}

		if rec.IsUnmapped() || rec.Cigar == "*" || rec.Seq == "*" || rec.Seq == "" {
			stats.noCall++
			if !polyaNoNA {
				emit(rec, chrom, "NA", "NA", "NA", "NA")
			}
			continue
		}

		mRNAPlus := !rec.IsReverse()
		if polyaAntisense {
			mRNAPlus = !mRNAPlus
		}
		strand := "+"
		if !mRNAPlus {
			strand = "-"
		}

		noCall := func() {
			stats.noCall++
			if !polyaNoNA {
				emit(rec, chrom, "NA", strand, "NA", "NA")
			}
		}

		if terminalHardClip(rec.Cigar, mRNAPlus) {
			noCall()
			continue
		}

		g := geomFor(len(rec.Seq), mRNAPlus)
		scanFrom := adapterScanStart(rec.Seq, g, aligner, polyaAdapter, polyaAdapterIdent, params.maxJunk)
		call := findPolyA(rec.Seq, g, params, scanFrom)
		if !call.found {
			noCall()
			continue
		}

		pos, src, ok := polyaRefPos(rec, call.start, g.step)
		if !ok {
			noCall()
			continue
		}
		// A tail at the very start of a contig, or an alignment ending exactly at
		// a contig's end, can push the site off the contig.
		if ln, known := refLens[rec.RefName]; known {
			switch {
			case pos < 1:
				pos = 1
				stats.clamped++
			case pos > ln:
				pos = ln
				stats.clamped++
			}
		} else if pos < 1 {
			pos = 1
			stats.clamped++
		}

		stats.called++
		emit(rec, chrom, strconv.Itoa(pos), strand, strconv.Itoa(call.length), src)
	}

	if err := w.Flush(); err != nil {
		closeFn()
		return fmt.Errorf("writing output: %w", err)
	}
	if err := closeFn(); err != nil {
		return err
	}

	fmt.Fprintf(os.Stderr, "Reads: %d, poly(A) called: %d, no call: %d\n",
		stats.total, stats.called, stats.noCall)
	if stats.skipped > 0 {
		fmt.Fprintf(os.Stderr, "Skipped %d secondary/supplementary alignments\n", stats.skipped)
	}
	if stats.clamped > 0 {
		fmt.Fprintf(os.Stderr, "Clamped %d position(s) to contig bounds\n", stats.clamped)
	}
	return nil
}

var ontPolyaCmd = &cobra.Command{
	GroupID:     "ontcmd",
	Annotations: map[string]string{"since": "v0.4.5"},
	Use:         "ont-polya <input.bam>",
	Short:       "Find poly(A)/cleavage sites from a strand-specific aligned BAM",
	Long: `For each read in a strand-specific aligned BAM, locate the mRNA 3' end, trace
back through the poly(A) tail to its first A base, and report that base's genomic
position -- the poly(A) / cleavage site.

The tail is detected with a sliding window: starting at the read's 3' terminus and
walking 5'-ward along the mRNA, bases are accepted while the trailing --window
bases keep an A fraction of at least --min-a-frac. The trace continues past the
soft-clip boundary into aligned bases, so genome-encoded A's that the aligner
absorbed into the alignment are counted as part of the tail.

Read orientation comes from FLAG 0x10, which is assumed to mean the read is sense
to the mRNA. Use --antisense for libraries where reads are antisense.

polya_pos is 1-based. Every read produces a row; NA marks a no-call. Secondary and
supplementary alignments are skipped: each read has exactly one primary, and
supplementary alignments are usually hard-clipped, so their tail is not present.

--polya-src reports where the tail's first base fell:
  aligned  mapped through the CIGAR; the aligner absorbed genome-encoded A's
  clip     in a soft clip, so the site is bounded at the alignment boundary
           rather than exactly mapped
  ins      in an insertion, anchored on the nearest reference base (+/-1)

--polya-length reports the observed tail in read bases, including genome-encoded
A's absorbed into the alignment. This measures something different from dorado's
pt tag, which is a raw-signal estimate, so the two will disagree -- emit both with
--polya-length --tag pt to compare them.

Examples:
  cgkit ont-polya aligned.bam
  cgkit ont-polya --min-a-frac 0.85 --min-len 15 aligned.bam
  cgkit ont-polya --polya-length --polya-src --tag pt -o sites.tsv aligned.bam
  cgkit ont-polya --antisense --adapter AGATCGGAAGAGC aligned.bam
  cgkit ont-polya --no-na -o sites.tsv.gz -t 4 aligned.bam`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if polyaWindow < 1 {
			return fmt.Errorf("--window must be at least 1")
		}
		if polyaMinLen < 1 {
			return fmt.Errorf("--min-len must be at least 1")
		}
		if polyaMinAFrac <= 0 || polyaMinAFrac > 1 {
			return fmt.Errorf("--min-a-frac must be greater than 0 and at most 1")
		}
		if polyaMaxJunk < 0 {
			return fmt.Errorf("--max-3p-junk cannot be negative")
		}
		if polyaAdapterIdent <= 0 || polyaAdapterIdent > 1 {
			return fmt.Errorf("--adapter-min-ident must be greater than 0 and at most 1")
		}
		return runPolya(args[0])
	},
}

var polyaOutput string
var polyaAntisense bool
var polyaMinAFrac float64
var polyaMinLen int
var polyaWindow int
var polyaMaxJunk int
var polyaAdapter string
var polyaAdapterIdent float64
var polyaShowLength bool
var polyaShowSource bool
var polyaTags []string
var polyaNoNA bool
var polyaThreads int
var polyaCramRef string

func init() {
	ontPolyaCmd.Flags().StringVarP(&polyaOutput, "output", "o", "-", "Output file path (default: stdout)")
	ontPolyaCmd.Flags().BoolVar(&polyaAntisense, "antisense", false, "Reads are antisense to the mRNA (invert the strand implied by FLAG 0x10)")
	ontPolyaCmd.Flags().Float64Var(&polyaMinAFrac, "min-a-frac", 0.8, "Minimum fraction of A bases within the sliding window")
	ontPolyaCmd.Flags().IntVar(&polyaMinLen, "min-len", 10, "Minimum poly(A) length (bases) required to make a call")
	ontPolyaCmd.Flags().IntVar(&polyaWindow, "window", 10, "Sliding window size (bases) for the A-fraction test")
	ontPolyaCmd.Flags().IntVar(&polyaMaxJunk, "max-3p-junk", 20, "Maximum non-A bases tolerated between the read 3' end and the tail")
	ontPolyaCmd.Flags().StringVar(&polyaAdapter, "adapter", "", "3' adapter sequence; when matched, the tail scan starts inside it")
	ontPolyaCmd.Flags().Float64Var(&polyaAdapterIdent, "adapter-min-ident", 0.75, "Minimum identity (matches / adapter length) for an adapter hit")
	ontPolyaCmd.Flags().BoolVar(&polyaShowLength, "polya-length", false, "Add a polya_len column (traced tail length, in read bases)")
	ontPolyaCmd.Flags().BoolVar(&polyaShowSource, "polya-src", false, "Add a polya_source column (aligned, clip, or ins)")
	ontPolyaCmd.Flags().Var(&tagArrayValue{values: &polyaTags}, "tag", "Append a column with this SAM tag's value; repeatable (e.g. --tag pt)")
	ontPolyaCmd.Flags().BoolVar(&polyaNoNA, "no-na", false, "Omit rows for reads with no poly(A) call")
	ontPolyaCmd.Flags().IntVarP(&polyaThreads, "threads", "t", 1, "Number of BGZF decompression threads")
	ontPolyaCmd.Flags().StringVar(&polyaCramRef, "cram-ref", "", "Reference FASTA for CRAM files")
}
