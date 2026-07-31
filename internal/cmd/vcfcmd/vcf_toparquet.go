package vcfcmd

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/compgenlab/cghts/varstore"
	"github.com/compgenlab/cghts/vcf"
	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/spf13/cobra"
)

var (
	vcfToParquetOut          string
	vcfToParquetRegion       string
	vcfToParquetMinDP        int
	vcfToParquetNoCallable   bool
	vcfToParquetPassing      bool
	vcfToParquetCompression  string
	vcfToParquetRowGroupSize int
	vcfToParquetVerbose      bool
	vcfToParquetForce        bool
)

var vcfToParquetCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": "v0.5.0"},
	Use:         "vcf-toparquet <input.vcf> [input2.vcf ...]",
	Short:       "Convert a VCF to a sparse Parquet genotype store",
	Long: `Convert a VCF into a columnar genotype store that keeps only the
alternate-allele calls, along with enough context to still tell a
confidently-called reference apart from a position that was never assayed.

Several inputs may be given, which is how whole-genome callsets usually ship --
one VCF per chromosome. They must carry exactly the same samples; differing
column order is remapped, since genotype columns are positional and getting that
wrong would silently attribute every genotype to the wrong person. A sample-set
mismatch is an error naming what differs.

Inputs must not overlap: a chromosome cannot be revisited once left, and
positions cannot go backwards within one. Overlapping inputs would write the
same site twice and split its AC/AN across two rows, so this is refused.

Give them in coordinate order for best query performance. Correctness does not
depend on it -- the answers are identical either way -- but the per-row-group
position statistics stay tighter, and a locus lookup then skips more of the
file. Measured on a two-chromosome store, supplying them out of order cost
about 1.8x on a locus query (166ms against 298ms).

Three files are written from --out, and they form one inseparable set:

  BASE.calls.parquet     one row per ALT-carrying genotype
  BASE.sites.parquet     one row per interrogated site, with AC/AN and counts
  BASE.regions.parquet   contiguous runs of adequately-covered sites, per sample

Ending --out with a "/" instead names a directory, which is created if needed,
and the members go inside it under their bare names:

  --out cohort/   ->  cohort/calls.parquet, cohort/sites.parquet, cohort/regions.parquet

That keeps the set as a single thing to copy, move or delete -- worth having,
since the three files are only meaningful together. vcf-varquery accepts either
form, and either member path within it.

Conversion refuses to overwrite an existing store: if any of the three members
is already present under --out, or if a prefix-form base names an existing
directory, it stops and asks for --force. Writing truncates all three, and a
half-replaced set is worse than either keeping or replacing the old one.

The sites file carries both allele counts (AC, AN) and sample counts
(n_carriers, n_called, n_lowdp). They are not interchangeable: a 1/1 genotype is
one carrier but two alt alleles, so AC >= n_carriers wherever a homozygote
occurs, and AN counts alleles without regard to depth while n_called counts
samples that cleared --min-dp. AF is exactly AC/AN. Both are computed over the
samples in this store rather than copied from the source's INFO fields, which
would be wrong after splitting a multiallelic record or converting a subset of a
cohort.

The sites file is not redundant with the calls. Deriving the site list from the
distinct loci in the calls only works when the store holds an entire joint
callset; over a subset of samples, every site where nobody in that subset
carries an ALT disappears, and a later query would report those positions as
never interrogated rather than as observed and reference.

Records are normalized to one variant per row: a multiallelic record is split
so each ALT allele gets its own rows. Within a split row the focal allele is
recoded to 1, reference stays 0, and any other alternate allele becomes "." --
so a 1/2 sample is correctly a carrier of both alleles. AD is taken per allele
(ad_ref is AD[0], ad_alt is that allele's own depth) rather than summed, since
depth supporting one alternate says nothing about another. Indels are NOT
left-aligned; normalize beforehand if the source is not already.

The regions file records, per sample, runs of catalog sites at which that
sample was successfully called at DP >= --min-dp. The interval form is only a
compression of that per-site fact; it makes NO claim about the bases between
those sites.

This bounds what the store can answer. A plain VCF reports variants and says
nothing whatsoever about any other position -- an unreported base was not
observed to be reference, it was simply never reported. The sites catalog is
therefore the exact boundary of what is knowable, and a query for a locus
outside it returns not-assayed for every sample rather than a set of reference
calls, even where run intervals appear to bracket it. Only a gVCF, whose
reference blocks carry END and MIN_DP, makes positive statements about spans and
could answer off-catalog positions.

  --out BASE            base name for the three output files, or DIR/ (required)
  --force               overwrite an existing store at --out
  --min-dp N            depth at or above which a site counts as callable
  --no-callable         proceed when the input has no DP field at all
  --passing             skip filtered records
  --compression C       zstd (default), snappy, or none
  --row-group-size N    rows per parquet row group
  -v, --verbose         report progress and a conversion summary on stderr,
                        including which of DP/GQ/AD the input actually carries`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if vcfToParquetOut == "" {
			return fmt.Errorf("you must specify a base output name with --out")
		}
		if vcfToParquetMinDP < 0 {
			return fmt.Errorf("--min-dp must not be negative")
		}
		if vcfToParquetRowGroupSize <= 0 {
			return fmt.Errorf("--row-group-size must be a positive number")
		}
		codec, err := varstore.CodecFor(vcfToParquetCompression)
		if err != nil {
			return err
		}

		// The first input fixes the sample roster; every later one must carry
		// the same people, though not necessarily in the same column order.
		first, err := openRecordSource(cmd, args[0], vcfToParquetRegion)
		if err != nil {
			return err
		}
		samples := first.header.Samples()
		first.close()
		if len(samples) == 0 {
			return fmt.Errorf("%s has no samples; a genotype store needs per-sample calls", args[0])
		}

		// Refuse to clobber an existing store before opening anything: the
		// writer truncates all three members, so this is the last moment the
		// previous one still exists.
		if err := varstore.CheckStoreTarget(vcfToParquetOut, vcfToParquetForce); err != nil {
			return err
		}
		// A base ending in "/" names a directory to put the members in.
		if err := varstore.EnsureStoreDir(vcfToParquetOut); err != nil {
			return err
		}

		// Before the writer exists, since it needs them at construction.
		contigs, err := collectContigs(cmd, args, vcfToParquetRegion)
		if err != nil {
			return err
		}
		if len(contigs) == 0 && vcfToParquetVerbose {
			fmt.Fprintf(cmd.ErrOrStderr(),
				"note: no ##contig lines in the input; a VCF exported from this store "+
					"cannot declare its reference\n")
		}

		w, err := varstore.NewWriter(vcfToParquetOut, varstore.WriterOpts{
			Codec:        codec,
			RowGroupSize: int64(vcfToParquetRowGroupSize),
			Samples:      samples,
			MinDP:        int32(vcfToParquetMinDP),
			NoCallable:   vcfToParquetNoCallable,
			Program:      buildinfo.String(),
			Command:      buildinfo.CommandLine(),
			Source:       strings.Join(args, ", "),
			Contigs:      contigs,
		})
		if err != nil {
			return err
		}

		started := time.Now()
		conv := &parquetConverter{
			w:          w,
			samples:    samples,
			minDP:      int32(vcfToParquetMinDP),
			noCallable: vcfToParquetNoCallable,
			runs:       make([]*callableRun, len(samples)),
			verbose:    vcfToParquetVerbose,
			progress:   cmd.ErrOrStderr(),
		}

		for _, path := range args {
			if err := convertOne(cmd, conv, path, samples); err != nil {
				w.Discard()
				return err
			}
		}
		if err := conv.finish(); err != nil {
			w.Discard()
			return err
		}
		if err := w.Close(); err != nil {
			return err
		}

		if conv.sawDP == 0 && !vcfToParquetNoCallable {
			return fmt.Errorf("no DP field found in %s, so callable regions cannot be built\n"+
				"       re-run with --no-callable to accept a store that cannot distinguish\n"+
				"       non-carrier from not-assayed", strings.Join(args, ", "))
		}
		if vcfToParquetVerbose {
			conv.report(cmd.ErrOrStderr(), vcfToParquetOut, time.Since(started))
		}
		fmt.Fprintf(cmd.ErrOrStderr(),
			"wrote %s: %d calls, %d sites, %d callable runs over %d samples\n",
			vcfToParquetOut, w.NCalls, w.NSites, w.NRegions, len(samples))
		return nil
	},
}

// samplePermutation maps this file's genotype columns onto the canonical
// sample order, and fails if the file does not carry exactly the same people.
//
// Genotype columns are addressed positionally, so getting this wrong does not
// error -- it silently attributes every genotype to the wrong person and
// produces entirely plausible output. Reordering is therefore remapped rather
// than merely tolerated: a bcftools merge or -S reorder is easy to do by
// accident, and remapping turns a silent corruption into a correct result.
func samplePermutation(canonical, got []string, path string) ([]int, bool, error) {
	if len(got) != len(canonical) {
		return nil, false, sampleMismatch(canonical, got, path)
	}
	index := make(map[string]int, len(canonical))
	for i, s := range canonical {
		index[s] = i
	}
	perm := make([]int, len(got))
	reordered := false
	for i, s := range got {
		j, ok := index[s]
		if !ok {
			return nil, false, sampleMismatch(canonical, got, path)
		}
		perm[i] = j
		if i != j {
			reordered = true
		}
	}
	return perm, reordered, nil
}

// sampleMismatch reports which samples differ, since "sample lists differ" on
// a 3,000-sample cohort is not an actionable message.
func sampleMismatch(canonical, got []string, path string) error {
	have := map[string]bool{}
	for _, s := range got {
		have[s] = true
	}
	want := map[string]bool{}
	for _, s := range canonical {
		want[s] = true
	}
	var missing, extra []string
	for _, s := range canonical {
		if !have[s] {
			missing = append(missing, s)
		}
	}
	for _, s := range got {
		if !want[s] {
			extra = append(extra, s)
		}
	}
	msg := fmt.Sprintf("%s does not carry the same samples as the first input (%d vs %d)",
		path, len(got), len(canonical))
	if len(missing) > 0 {
		msg += "\n       missing: " + summariseNames(missing)
	}
	if len(extra) > 0 {
		msg += "\n       unexpected: " + summariseNames(extra)
	}
	return fmt.Errorf("%s", msg)
}

// summariseNames lists a few names rather than thousands.
func summariseNames(names []string) string {
	const max = 6
	if len(names) <= max {
		return strings.Join(names, ", ")
	}
	return fmt.Sprintf("%s and %d more", strings.Join(names[:max], ", "), len(names)-max)
}

// convertOne streams one input into the store.
func convertOne(cmd *cobra.Command, conv *parquetConverter, path string, canonical []string) error {
	src, err := openRecordSource(cmd, path, vcfToParquetRegion)
	if err != nil {
		return err
	}
	defer src.close()

	perm, reordered, err := samplePermutation(canonical, src.header.Samples(), path)
	if err != nil {
		return err
	}
	conv.perm = perm

	if conv.verbose {
		note := ""
		if reordered {
			note = "  (sample columns reordered to match the first input)"
		}
		fmt.Fprintf(conv.progress, "reading %s (%d samples)%s\n", path, len(canonical), note)
	}
	conv.nFiles++

	for {
		rec, err := src.next()
		if err == io.EOF {
			return nil
		}
		if err != nil {
			return err
		}
		if vcfToParquetPassing && rec.IsFiltered() {
			conv.nFiltered++
			continue
		}
		if err := conv.checkOrder(rec, path); err != nil {
			return err
		}
		if err := conv.record(rec); err != nil {
			return err
		}
	}
}

// callableRun is an in-progress run of covered sites for one sample.
type callableRun struct {
	start  int32
	last   int32
	nSites int32
}

// parquetConverter turns a stream of VCF records into store rows, holding only
// one open run per sample so memory does not grow with the input.
type parquetConverter struct {
	w          *varstore.Writer
	samples    []string
	minDP      int32
	noCallable bool
	runs       []*callableRun
	curChrom   string
	sawDP      int64

	// perm maps the current file's genotype columns onto canonical sample
	// indices; nil means the columns are already in canonical order.
	perm []int

	// ordering cursor, to keep the concatenation coordinate sorted
	lastPos   int32
	seenChrom map[string]bool
	nFiles    int

	// verbose reporting
	verbose       bool
	progress      io.Writer
	nRecords      int64
	nFiltered     int64
	nMultiAllelic int64
	nExtraRows    int64
	sawGQ         int64
	sawAD         int64
	nGenotypes    int64
	nBelowDP      int64
	nNoCall       int64
	chroms        []string
}

// note records a per-chromosome transition for the verbose summary.
func (c *parquetConverter) note(chrom string) {
	c.chroms = append(c.chroms, chrom)
	if c.verbose {
		fmt.Fprintf(c.progress, "  %s: starting\n", chrom)
	}
}

// tick emits periodic progress, which matters because a whole-chromosome
// conversion streams for minutes with nothing else to show for it.
func (c *parquetConverter) tick() {
	const every = 100_000
	if c.verbose && c.nRecords%every == 0 {
		fmt.Fprintf(c.progress, "  %s: %d records, %d calls so far\n",
			c.curChrom, c.nRecords, c.w.NCalls)
	}
}

// record splits one VCF record into per-allele calls, a catalog entry per
// allele, and callable-run bookkeeping.
func (c *parquetConverter) record(rec *vcf.VcfRecord) error {
	chrom := rec.Chrom
	if chrom != c.curChrom {
		// Runs cannot span chromosomes.
		if err := c.closeRuns(); err != nil {
			return err
		}
		c.curChrom = chrom
		c.note(chrom)
	}
	c.nRecords++
	c.tick()

	alts := rec.Alt()
	pos := int32(rec.Pos)
	nAlts := len(alts)
	if nAlts > 1 {
		c.nMultiAllelic++
		c.nExtraRows += int64(nAlts - 1)
	}
	carriers := make([]int32, nAlts)
	acCounts := make([]int32, nAlts)
	var an int32
	var nLowDP, nCalled int32

	n := rec.NumSamples()
	if n > len(c.samples) {
		n = len(c.samples)
	}
	for i := 0; i < n; i++ {
		sf, err := varstore.ReadSample(rec, i)
		if err != nil {
			return fmt.Errorf("%w (%s:%d)", err, rec.Chrom, rec.Pos)
		}
		c.nGenotypes++
		if sf.DP != varstore.Missing {
			c.sawDP++
		}
		if sf.GQ != varstore.Missing {
			c.sawGQ++
		}
		if sf.AD != "" {
			c.sawAD++
		}
		if !varstore.HasCall(sf.GT) {
			c.nNoCall++
		} else if sf.DP != varstore.Missing && sf.DP < c.minDP {
			c.nBelowDP++
		}

		// Allele counts come straight from GT and are deliberately outside the
		// --no-callable guard below: AC/AN are properties of the genotypes, not
		// of coverage, so they stay meaningful even for a source with no DP.
		an += varstore.AddAlleleCounts(sf.GT, acCounts)

		// Coverage bookkeeping is per site, not per allele. A site counts as
		// callable only when the caller actually made a call there AND depth
		// clears the threshold; "./." at high depth is a declined call, not a
		// covered one.
		if !c.noCallable {
			si := c.sampleAt(i)
			if varstore.HasCall(sf.GT) && sf.DP != varstore.Missing && sf.DP >= c.minDP {
				nCalled++
				if r := c.runs[si]; r != nil {
					r.last = pos
					r.nSites++
				} else {
					c.runs[si] = &callableRun{start: pos, last: pos, nSites: 1}
				}
			} else {
				nLowDP++
				if err := c.emitRun(si); err != nil {
					return err
				}
			}
		}

		name := c.samples[c.sampleAt(i)]
		for j, alt := range alts {
			call, ok := varstore.CallFor(rec, name, sf, j+1, alt)
			if !ok {
				continue
			}
			carriers[j]++
			if err := c.w.WriteCall(call); err != nil {
				return err
			}
		}
	}

	for j, alt := range alts {
		if err := c.w.WriteSite(varstore.Site{
			Chrom:     chrom,
			Pos:       pos,
			Ref:       rec.Ref,
			Alt:       alt,
			AC:        acCounts[j],
			AN:        an,
			NCarriers: carriers[j],
			NLowDP:    nLowDP,
			NCalled:   nCalled,
		}); err != nil {
			return err
		}
	}
	return nil
}

// emitRun writes and clears sample i's open run, if any.
func (c *parquetConverter) emitRun(i int) error {
	r := c.runs[i]
	if r == nil {
		return nil
	}
	c.runs[i] = nil
	return c.w.WriteRegion(varstore.CalledSiteRun{
		SampleID: c.samples[i],
		Chrom:    c.curChrom,
		Start:    r.start,
		End:      r.last,
		NSites:   r.nSites,
	})
}

// closeRuns flushes every open run, at a chromosome change or end of input.
func (c *parquetConverter) closeRuns() error {
	for i := range c.runs {
		if err := c.emitRun(i); err != nil {
			return err
		}
	}
	return nil
}

// finish flushes any runs still open at end of input.
func (c *parquetConverter) finish() error { return c.closeRuns() }

// report writes the verbose conversion summary.
//
// The field-presence section is the part worth having. A gate can only act on a
// field the data carries, and --min-gq over GQ-less input admits everything
// rather than rejecting it. That is deliberate -- absent quality is not evidence
// of poor quality -- but it means a filter can silently do nothing, so a store
// built from such input should say so at the point it is created.
func (c *parquetConverter) report(out io.Writer, base string, elapsed time.Duration) {
	pct := func(n int64) string {
		if c.nGenotypes == 0 {
			return "0%"
		}
		return fmt.Sprintf("%.1f%%", 100*float64(n)/float64(c.nGenotypes))
	}

	fmt.Fprintf(out, "\ninput\n")
	fmt.Fprintf(out, "  records read          %d\n", c.nRecords)
	if c.nFiltered > 0 {
		fmt.Fprintf(out, "  skipped (--passing)   %d\n", c.nFiltered)
	}
	fmt.Fprintf(out, "  multiallelic records  %d (split into %d extra rows)\n",
		c.nMultiAllelic, c.nExtraRows)
	fmt.Fprintf(out, "  chromosomes           %s\n", strings.Join(c.chroms, ", "))
	fmt.Fprintf(out, "  genotypes examined    %d\n", c.nGenotypes)

	fmt.Fprintf(out, "\nfields present (a gate can only act on a field the data has)\n")
	for _, f := range []struct {
		name string
		n    int64
		note string
	}{
		{"DP", c.sawDP, "callable runs and --min-dp depend on this"},
		{"GQ", c.sawGQ, "--min-gq depends on this"},
		{"AD", c.sawAD, "per-allele depths"},
	} {
		state := pct(f.n)
		if f.n == 0 {
			state = "ABSENT -- " + f.note + " will have no effect"
		}
		fmt.Fprintf(out, "  %-3s %s\n", f.name, state)
	}

	fmt.Fprintf(out, "\ncoverage at --min-dp %d\n", c.minDP)
	if c.noCallable {
		fmt.Fprintf(out, "  not tracked (--no-callable)\n")
	} else {
		fmt.Fprintf(out, "  no call made          %d  (%s)\n", c.nNoCall, pct(c.nNoCall))
		fmt.Fprintf(out, "  called but below DP   %d  (%s)\n", c.nBelowDP, pct(c.nBelowDP))
	}

	fmt.Fprintf(out, "\noutput\n")
	for _, f := range []struct {
		path string
		n    int64
	}{
		{varstore.CallsPath(base), c.w.NCalls},
		{varstore.SitesPath(base), c.w.NSites},
		{varstore.RegionsPath(base), c.w.NRegions},
	} {
		size := int64(-1)
		if st, err := os.Stat(f.path); err == nil {
			size = st.Size()
		}
		fmt.Fprintf(out, "  %-24s %9d rows  %10d bytes\n", filepath.Base(f.path), f.n, size)
	}
	fmt.Fprintf(out, "  elapsed               %s\n", elapsed.Round(time.Millisecond))
}

func init() {
	f := vcfToParquetCmd.Flags()
	f.StringVar(&vcfToParquetOut, "out", "", "Base output name; BASE.calls.parquet etc, or DIR/ for DIR/calls.parquet etc (the directory is created)")
	f.StringVar(&vcfToParquetRegion, "region", "", "Only variants in this 1-based region (chrom:start-end, or chrom); requires a tabix-indexed file")
	f.IntVar(&vcfToParquetMinDP, "min-dp", 10, "Minimum DP for a site to count as callable for a sample")
	f.BoolVar(&vcfToParquetNoCallable, "no-callable", false, "Accept a source with no DP field; callable regions will be empty")
	f.BoolVar(&vcfToParquetPassing, "passing", false, "Only convert passing variants")
	f.StringVar(&vcfToParquetCompression, "compression", "zstd", "Parquet compression: zstd, snappy, or none")
	f.IntVar(&vcfToParquetRowGroupSize, "row-group-size", 250000, "Rows per parquet row group")
	f.BoolVarP(&vcfToParquetVerbose, "verbose", "v", false, "Report progress and a conversion summary on stderr")
	f.BoolVar(&vcfToParquetForce, "force", false, "Overwrite an existing store at --out")
}

// sampleAt maps a genotype column of the file being read to its canonical
// sample index.
func (c *parquetConverter) sampleAt(col int) int {
	if c.perm == nil || col >= len(c.perm) {
		return col
	}
	return c.perm[col]
}

// checkOrder enforces that the inputs, concatenated, stay coordinate sorted.
//
// This is one rule serving three purposes: it keeps parquet's per-row-group
// min/max on pos tight, which is what makes locus lookups prune; it catches
// inputs supplied in the wrong order; and it rejects overlapping inputs, which
// would otherwise write duplicate sites and split AC/AN across two rows for the
// same variant.
func (c *parquetConverter) checkOrder(rec *vcf.VcfRecord, path string) error {
	if c.seenChrom == nil {
		c.seenChrom = map[string]bool{}
	}
	chrom, pos := rec.Chrom, int32(rec.Pos)
	if chrom == c.curChrom {
		if pos < c.lastPos {
			return fmt.Errorf("%s is not coordinate sorted at %s:%d (previous record was %d)\n"+
				"       inputs must be sorted, and must not overlap each other",
				path, chrom, pos, c.lastPos)
		}
		c.lastPos = pos
		return nil
	}
	if c.seenChrom[chrom] {
		return fmt.Errorf("%s returns to %s after another chromosome\n"+
			"       inputs must be in coordinate order and must not overlap; "+
			"pass them one chromosome at a time, in order", path, chrom)
	}
	c.seenChrom[chrom] = true
	c.lastPos = pos
	return nil
}

// collectContigs gathers the ##contig lines of every input, so the store can be
// exported back to VCF with a header that says which reference it came from.
//
// The UNION, not the first input's. A whole-genome callset usually ships one VCF
// per chromosome, and such a file often declares only its own contig -- taking the
// first input's list would silently lose every other chromosome, which is exactly
// the case multi-input conversion exists to serve.
//
// Lines are kept verbatim, so lengths and any extra fields survive rather than
// being reconstructed approximately.
func collectContigs(cmd *cobra.Command, inputs []string, region string) ([]string, error) {
	type origin struct {
		line string
		from string
	}
	byID := map[string]origin{}
	var order []string

	for _, path := range inputs {
		src, err := openRecordSource(cmd, path, region)
		if err != nil {
			return nil, err
		}
		names := src.header.ContigNames()
		for _, id := range names {
			def, ok := src.header.ContigDef(id)
			if !ok {
				continue
			}
			line := def.String()
			prev, seen := byID[id]
			if !seen {
				byID[id] = origin{line: line, from: path}
				order = append(order, id)
				continue
			}
			// Same contig declared differently by two inputs. A differing length
			// means the inputs were called against different references, which would
			// make one store out of two incompatible callsets -- refused for the same
			// reason a differing sample set is.
			if prev.line != line {
				src.close()
				return nil, fmt.Errorf("inputs disagree about contig %s:\n"+
					"       %s: %s\n"+
					"       %s: %s\n"+
					"       a differing length means these were called against different "+
					"references", id, prev.from, prev.line, path, line)
			}
		}
		src.close()
	}

	out := make([]string, 0, len(order))
	for _, id := range order {
		out = append(out, byID[id].line)
	}
	return out, nil
}
