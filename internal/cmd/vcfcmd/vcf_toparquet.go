package vcfcmd

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/compgenlab/cghts/vcf"
	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/compgenlab/cgkit/internal/varstore"
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
	Use:         "vcf-toparquet <input.vcf>",
	Short:       "Convert a VCF to a sparse Parquet genotype store",
	Long: `Convert a VCF into a columnar genotype store that keeps only the
alternate-allele calls, along with enough context to still tell a
confidently-called reference apart from a position that was never assayed.

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

		src, err := openRecordSource(cmd, args[0], vcfToParquetRegion)
		if err != nil {
			return err
		}
		defer src.close()

		samples := src.header.Samples()
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

		w, err := varstore.NewWriter(vcfToParquetOut, varstore.WriterOpts{
			Codec:        codec,
			RowGroupSize: int64(vcfToParquetRowGroupSize),
			Samples:      samples,
			MinDP:        int32(vcfToParquetMinDP),
			NoCallable:   vcfToParquetNoCallable,
			Program:      buildinfo.String(),
			Command:      buildinfo.CommandLine(),
			Source:       args[0],
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
		if vcfToParquetVerbose {
			fmt.Fprintf(cmd.ErrOrStderr(), "reading %s (%d samples)\n", args[0], len(samples))
		}

		for {
			rec, err := src.next()
			if err == io.EOF {
				break
			}
			if err != nil {
				w.Close()
				return err
			}
			if vcfToParquetPassing && rec.IsFiltered() {
				conv.nFiltered++
				continue
			}
			if err := conv.record(rec); err != nil {
				w.Close()
				return err
			}
		}
		if err := conv.finish(); err != nil {
			w.Close()
			return err
		}
		if err := w.Close(); err != nil {
			return err
		}

		if conv.sawDP == 0 && !vcfToParquetNoCallable {
			return fmt.Errorf("no DP field found in %s, so callable regions cannot be built\n"+
				"       re-run with --no-callable to accept a store that cannot distinguish\n"+
				"       non-carrier from not-assayed", args[0])
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
			if varstore.HasCall(sf.GT) && sf.DP != varstore.Missing && sf.DP >= c.minDP {
				nCalled++
				if r := c.runs[i]; r != nil {
					r.last = pos
					r.nSites++
				} else {
					c.runs[i] = &callableRun{start: pos, last: pos, nSites: 1}
				}
			} else {
				nLowDP++
				if err := c.emitRun(i); err != nil {
					return err
				}
			}
		}

		for j, alt := range alts {
			call, ok := varstore.CallFor(rec, c.samples[i], sf, j+1, alt)
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
