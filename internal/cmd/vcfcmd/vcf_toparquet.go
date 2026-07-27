package vcfcmd

import (
	"fmt"
	"io"

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
)

var vcfToParquetCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": "v0.5.0"},
	Use:         "vcf-toparquet <input.vcf>",
	Short:       "Convert a VCF to a sparse Parquet genotype store",
	Long: `Convert a VCF into a columnar genotype store that keeps only the
alternate-allele calls, along with enough context to still tell a
confidently-called reference apart from a position that was never assayed.

Three files are written from --out BASE, and they form one inseparable set:

  BASE.calls.parquet     one row per ALT-carrying genotype
  BASE.sites.parquet     one row per interrogated site, sample-independent
  BASE.regions.parquet   contiguous runs of adequately-covered sites, per sample

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

Callable regions are runs of consecutive variant sites at which a sample had
DP >= --min-dp. The span between two in-run sites is assumed callable, not
observed: a plain VCF says nothing between its records. Only a gVCF, whose
reference blocks carry END and MIN_DP, would make these true observations.

  --out BASE            base name for the three output files (required)
  --min-dp N            depth at or above which a site counts as callable
  --no-callable         proceed when the input has no DP field at all
  --passing             skip filtered records
  --compression C       zstd (default), snappy, or none
  --row-group-size N    rows per parquet row group`,
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

		conv := &parquetConverter{
			w:          w,
			samples:    samples,
			minDP:      int32(vcfToParquetMinDP),
			noCallable: vcfToParquetNoCallable,
			runs:       make([]*callableRun, len(samples)),
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
}

// record splits one VCF record into per-allele calls, a catalog entry per
// allele, and callable-run bookkeeping.
func (c *parquetConverter) record(rec *vcf.VcfRecord) error {
	chrom := varstore.NormChrom(rec.Chrom)
	if chrom != c.curChrom {
		// Runs cannot span chromosomes.
		if err := c.closeRuns(); err != nil {
			return err
		}
		c.curChrom = chrom
	}

	alts := rec.Alt()
	pos := int32(rec.Pos)
	nAlts := len(alts)
	carriers := make([]int32, nAlts)
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
		if sf.DP != varstore.Missing {
			c.sawDP++
		}

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
	return c.w.WriteRegion(varstore.CallableRegion{
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

func init() {
	f := vcfToParquetCmd.Flags()
	f.StringVar(&vcfToParquetOut, "out", "", "Base output name (outputs are BASE.calls.parquet, BASE.sites.parquet, BASE.regions.parquet)")
	f.StringVar(&vcfToParquetRegion, "region", "", "Only variants in this 1-based region (chrom:start-end, or chrom); requires a tabix-indexed file")
	f.IntVar(&vcfToParquetMinDP, "min-dp", 10, "Minimum DP for a site to count as callable for a sample")
	f.BoolVar(&vcfToParquetNoCallable, "no-callable", false, "Accept a source with no DP field; callable regions will be empty")
	f.BoolVar(&vcfToParquetPassing, "passing", false, "Only convert passing variants")
	f.StringVar(&vcfToParquetCompression, "compression", "zstd", "Parquet compression: zstd, snappy, or none")
	f.IntVar(&vcfToParquetRowGroupSize, "row-group-size", 250000, "Rows per parquet row group")
}
