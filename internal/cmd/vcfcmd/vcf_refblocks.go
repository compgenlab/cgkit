package vcfcmd

import (
	"errors"
	"fmt"
	"io"
	"slices"

	"github.com/compgenlab/cghts/varstore"
	_ "github.com/compgenlab/cghts/varstore/sinks3" // --out s3://...
	"github.com/compgenlab/cghts/vcf"
	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/spf13/cobra"
)

var (
	vcfRefBlocksOut         string
	vcfRefBlocksRoster      string
	vcfRefBlocksMinDP       int
	vcfRefBlocksDepthRatio  float64
	vcfRefBlocksBands       []int
	vcfRefBlocksMaxGap      int
	vcfRefBlocksPassthrough bool
	vcfRefBlocksForce       bool
)

var vcfRefBlocksCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": "v0.6.11"},
	Use:         "vcf-refblocks <sample.g.vcf.gz>",
	Short:       "Extract one sample's reference blocks from a gVCF, for a varstore",
	Long: `Read one sample's gVCF and write its reference blocks as a FRAGMENT: the
positive coverage claims a store's ref_blocks table is built from.

A store converted from a joint-called pVCF can only answer at its catalog. A
variant private to another batch, or found after conversion, is not-assayed for
everyone -- the pVCF never said who was reference there. A gVCF did, and this is
how that is captured. Genotypes still come from the pVCF; a fragment holds no
ALT calls, only spans where the sample was observed to be reference.

Run it once per sample, independently -- it is single-sample and streaming, and
the fragments are merged into a store afterwards (vcf-varstore attach-refblocks).

  cgkit vcf-refblocks S0001.g.vcf.gz --roster cohort/ --out frags/S0001

What is kept:

  reference blocks   records whose only alternates are <NON_REF> / <*>, with
                     END and MIN_DP. Depth is MIN_DP (falling back to DP); GQ
                     (RGQ) is stored but never gated on.
  below --min-dp     dropped -- uncovered at the gate. For exome data plus
                     low-pass WGS this leaves essentially the exome.

What is never merged across: a VARIANT RECORD -- any record with a real
alternate, filtered or not -- and a block that declined to call (./.). Those
are positions where the sample was not, or may not have been, reference; their
gap is what keeps a block from claiming reference where the sample carries
something, so no tolerance bridges them.

Coarsening, because a GATK gVCF breaks a block at every GQ change and a store
never gates on GQ:

  depth bands   adjacent blocks in the same depth band merge, keeping the
                lowest MIN_DP. The bands default to the roster store's own, so
                a query gate means the same thing for blocks and catalog runs.
  --max-gap N   merge across uncovered stretches of at most N bases (default
                10). This weakens the claim to "covered at min_dp, with no
                uncovered run longer than N", so N is recorded.
  --passthrough keep the source's blocks as written: no banding, no merging.
                For targeted or exome work, where fidelity is cheap.

--roster names the store the fragment is for. The gVCF's sample must be on its
roster, and --min-dp and the bands default to the store's; giving different
ones is refused, since a fragment that cannot be attached is wasted work.

The output is a directory with one Parquet file per chromosome and
fragment.json, written LAST, which marks it complete. A run that dies leaves no
marker, so it is safe to rerun; an existing finished fragment is kept unless
--force.

  --out DIR             the fragment directory (required); may be s3://...
  --roster STORE        the store the fragment will be attached to (required)
  --min-dp N            drop blocks below this depth (default: the store's)
  --depth-ratio R       depth bands at min-dp x R^k (default: the store's bands)
  --depth-bands D,...   explicit depth bands instead
  --max-gap N           longest uncovered stretch a merge may span (default 10)
  --passthrough         no banding and no merging
  --force               replace an existing finished fragment`,
	Args: cobra.ExactArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		ctx := cmd.Context()
		in := args[0]
		if vcfRefBlocksOut == "" {
			return fmt.Errorf("--out is required")
		}
		if vcfRefBlocksRoster == "" {
			return fmt.Errorf("--roster is required: a fragment is built for one store, " +
				"and its sample and settings are checked against that store before any work is done")
		}
		if vcfRefBlocksMaxGap < 0 {
			return fmt.Errorf("--max-gap %d is negative", vcfRefBlocksMaxGap)
		}
		if !varstore.CanWrite(vcfRefBlocksOut) {
			return fmt.Errorf("cannot write to %s: no transport for its scheme is linked into this build",
				vcfRefBlocksOut)
		}

		// The roster first, and the source's header second: both are cheap, and
		// both can refuse the job before the gVCF is read at all.
		roster, params, err := rosterOf(cmd, vcfRefBlocksRoster)
		if err != nil {
			return err
		}
		minDP, bands, ratio, err := refBlocksSettings(cmd, params)
		if err != nil {
			return err
		}

		r, err := vcf.OpenVcfFile(ctx, in)
		if err != nil {
			return err
		}
		defer r.Close()
		h, err := r.Header()
		if err != nil {
			return err
		}
		if n := len(h.Samples()); n != 1 {
			return fmt.Errorf("%s has %d samples; a gVCF fragment is single-sample by construction", in, n)
		}
		sample := h.Samples()[0]
		if !slices.Contains(roster, sample) {
			return fmt.Errorf("%s is not among the %d samples of %s; its blocks could never be "+
				"attached, so nothing was written", sample, len(roster), vcfRefBlocksRoster)
		}

		sink, err := varstore.OpenSink(vcfRefBlocksOut)
		if err != nil {
			return err
		}
		if _, done, err := sink.Stat(varstore.FragmentMetaFile); err != nil {
			return err
		} else if done && !vcfRefBlocksForce {
			return fmt.Errorf("%s already holds a finished fragment; use --force to replace it",
				vcfRefBlocksOut)
		}

		fw, err := varstore.NewFragmentWriter(vcfRefBlocksOut, varstore.FragmentMeta{
			Program:     buildinfo.String(),
			Command:     buildinfo.CommandLine(),
			Source:      in,
			SampleID:    sample,
			MinDP:       minDP,
			DepthBands:  bands,
			DepthRatio:  ratio,
			MaxGap:      int32(vcfRefBlocksMaxGap),
			Passthrough: vcfRefBlocksPassthrough,
		})
		if err != nil {
			return err
		}
		stats, err := extractRefBlocks(r, sample, fw, varstore.Coarsener{
			MinDP:       minDP,
			Bands:       bands,
			MaxGap:      int32(vcfRefBlocksMaxGap),
			Passthrough: vcfRefBlocksPassthrough,
		})
		if err == nil {
			err = fw.Finish()
		}
		if err != nil {
			if derr := fw.Discard(); derr != nil {
				return errors.Join(err, fmt.Errorf("and removing the partial fragment failed: %w", derr))
			}
			return err
		}
		m := fw.Meta()
		var rows int64
		for _, c := range m.Chroms {
			rows += c.Rows
		}
		fmt.Fprintf(cmd.ErrOrStderr(),
			"wrote %s: %s, %d source blocks -> %d spans on %d chromosome(s); %d below --min-dp %d, %d variant/no-call barriers\n",
			vcfRefBlocksOut, sample, stats.blocks, rows, len(m.Chroms), stats.dropped, minDP, stats.barriers)
		return nil
	},
}

// rosterOf reads a store's samples and conversion settings without opening its
// tables: the manifest is all a fragment needs to agree with.
func rosterOf(cmd *cobra.Command, locator string) ([]string, varstore.ManifestParams, error) {
	store, err := varstore.Open(cmd.Context(), locator, "")
	if err != nil {
		return nil, varstore.ManifestParams{}, err
	}
	defer store.Close()
	samples, err := store.Samples()
	if err != nil {
		return nil, varstore.ManifestParams{}, err
	}
	switch s := store.(type) {
	case *varstore.VarStore:
		return samples, s.StoreManifest().Params, nil
	case *varstore.ParquetVolume:
		return samples, s.VolumeManifest().Params, nil
	}
	return nil, varstore.ManifestParams{}, fmt.Errorf(
		"--roster %s is not a varstore; fragments are attached to a store, not a VCF", locator)
}

// refBlocksSettings resolves the gate and bands a fragment is built with.
//
// They default to the store's own, and a value given explicitly must still agree
// with it. A fragment banded or gated differently from its store would be
// refused at attach -- after the job that built it, one of perhaps ten
// thousand, has already run.
func refBlocksSettings(cmd *cobra.Command, p varstore.ManifestParams) (int32, []int32, float64, error) {
	minDP := p.MinDP
	if cmd.Flags().Changed("min-dp") {
		if int32(vcfRefBlocksMinDP) != p.MinDP {
			return 0, nil, 0, fmt.Errorf("--min-dp %d disagrees with the store's %d; "+
				"blocks gated differently from the catalog runs answer a gate differently", vcfRefBlocksMinDP, p.MinDP)
		}
	}
	bands, ratio := p.DepthBands, p.DepthRatio
	switch {
	case vcfRefBlocksPassthrough:
		return minDP, bands, ratio, nil // recorded, not applied
	case cmd.Flags().Changed("depth-bands") && cmd.Flags().Changed("depth-ratio"):
		return 0, nil, 0, fmt.Errorf("--depth-bands and --depth-ratio are alternatives; give one")
	case cmd.Flags().Changed("depth-bands"):
		bands = make([]int32, 0, len(vcfRefBlocksBands))
		for _, b := range vcfRefBlocksBands {
			bands = append(bands, int32(b))
		}
		ratio = 0
	case cmd.Flags().Changed("depth-ratio"):
		ratio = vcfRefBlocksDepthRatio
		if ratio == 0 {
			bands = nil
		} else {
			start := minDP
			if start < 1 {
				start = 1
			}
			var err error
			if bands, err = varstore.GeometricDepthBands(start, ratio, varstore.MaxDepthBand); err != nil {
				return 0, nil, 0, fmt.Errorf("--depth-ratio: %w", err)
			}
		}
	default:
		return minDP, bands, ratio, nil
	}
	if !slices.Equal(bands, p.DepthBands) {
		return 0, nil, 0, fmt.Errorf("depth bands %v disagree with the store's %v; "+
			"a block's band would mean something different from a catalog run's", bands, p.DepthBands)
	}
	return minDP, bands, ratio, nil
}

type refBlockStats struct {
	blocks, dropped, barriers int64
}

// extractRefBlocks streams the gVCF into the coarsener and the coarsener into
// the fragment.
func extractRefBlocks(r *vcf.VcfReader, sample string, fw *varstore.FragmentWriter, co varstore.Coarsener) (refBlockStats, error) {
	var st refBlockStats
	co.Emit = fw.Write
	for {
		rec, err := r.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			return st, err
		}
		if !rec.IsRefBlock() {
			// A variant record, real alternate present, filtered or not.
			st.barriers++
			if err := co.Barrier(); err != nil {
				return st, err
			}
			continue
		}
		sf, err := varstore.ReadSample(rec, 0)
		if err != nil {
			return st, err
		}
		if !varstore.IsHomRef(sf.GT) {
			// A block that declined to call says nothing about reference.
			st.barriers++
			if err := co.Barrier(); err != nil {
				return st, err
			}
			continue
		}
		st.blocks++
		// Unknown depth vouches for nothing, so it reads as zero and is dropped
		// by any positive gate.
		dp, _ := rec.BlockDepth(0)
		gq := varstore.Missing
		if v, ok := rec.BlockRGQ(0); ok {
			gq = v
		}
		if dp < co.MinDP {
			st.dropped++
		}
		if err := co.Block(varstore.RefBlock{
			SampleID: sample,
			Chrom:    rec.Chrom,
			Start:    int32(rec.Pos),
			End:      int32(rec.RefSpanEnd()), // 0-based exclusive == 1-based last base
			MinDP:    dp,
			GQ:       gq,
		}); err != nil {
			return st, err
		}
	}
	return st, co.Flush()
}

func init() {
	f := vcfRefBlocksCmd.Flags()
	f.StringVar(&vcfRefBlocksOut, "out", "", "Fragment directory, created if needed")
	f.StringVar(&vcfRefBlocksRoster, "roster", "", "The store the fragment will be attached to")
	f.IntVar(&vcfRefBlocksMinDP, "min-dp", 10, "Drop blocks below this depth (default: the store's)")
	f.Float64Var(&vcfRefBlocksDepthRatio, "depth-ratio", varstore.DefaultDepthRatio,
		"Depth bands at min-dp x ratio^k (default: the store's bands)")
	f.IntSliceVar(&vcfRefBlocksBands, "depth-bands", nil, "Explicit depth bands (default: the store's)")
	f.IntVar(&vcfRefBlocksMaxGap, "max-gap", 10, "Longest uncovered stretch, in bases, a merge may span")
	f.BoolVar(&vcfRefBlocksPassthrough, "passthrough", false, "Keep the source's blocks: no banding, no merging")
	f.BoolVar(&vcfRefBlocksForce, "force", false, "Replace an existing finished fragment")
}
