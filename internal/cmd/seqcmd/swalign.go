package seqcmd

import (
	"fmt"

	"github.com/compgenlab/cghts/align"
	"github.com/compgenlab/cghts/seqio"
	"github.com/spf13/cobra"
)

// swalignCmd implements the seq-pairwise command: Smith-Waterman pairwise alignment.
var swalignCmd = &cobra.Command{
	GroupID:     "seqcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "seq-pairwise query target",
	Short:       "Align the two given sequences",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) != 2 {
			cmd.Help()
			return nil
		}

		if swalignUseClipping && swalignGlobal {
			return fmt.Errorf("cannot use clipping penalties with global alignment")
		}

		// --gap-open/--gap-extend set both strands' penalties; the -ins/-del
		// forms override whichever they name. Both were registered and never
		// read, so passing --gap-open changed nothing at all.
		openIns, openDel := swalignGapOpenIns, swalignGapOpenDel
		extendIns, extendDel := swalignGapExtendIns, swalignGapExtendDel
		if f := cmd.Flags(); f.Changed("gap-open") {
			if !f.Changed("gap-open-ins") {
				openIns = swalignGapOpen
			}
			if !f.Changed("gap-open-del") {
				openDel = swalignGapOpen
			}
		}
		if f := cmd.Flags(); f.Changed("gap-extend") {
			if !f.Changed("gap-extend-ins") {
				extendIns = swalignGapExtend
			}
			if !f.Changed("gap-extend-del") {
				extendDel = swalignGapExtend
			}
		}

		opts := align.DnaAlignmentDefaults().
			ScoringMatrix(align.MatchMismatchScoring(swalignMatchScore, swalignMismatchPenalty)).
			GapPenaltyIns(openIns, extendIns).
			GapPenaltyDel(openDel, extendDel).
			HomopolymerDiscount(swalignHPOpenScale, swalignHPOpenCap, swalignHPExtendScale, swalignHPExtendCap).
			Verbose(swalignVerbose)

		if swalignUseClipping {
			opts = opts.ClippingPenalty(swalignClipOpen, swalignClipExtend)
		} else {
			opts = opts.ClippingDisable()
		}

		var sw align.PairwiseAligner
		if swalignGlobal {
			sw = align.NewGlobalAligner(opts)
		} else {
			sw = align.NewLocalAligner(opts)
		}
		query := seqio.NewStringSeq(args[0], "query")
		target := seqio.NewStringSeq(args[1], "target")
		aln1 := sw.Align(query.FullSeq(), target.FullSeq())
		aln2 := sw.Align(query.FullSeq().RevComp(), target.FullSeq())
		if aln1.Score >= aln2.Score {
			fmt.Println(aln1.String())
		} else {
			fmt.Println(aln2.String())
		}
		return nil
	},
}

var (
	swalignMatchScore      int
	swalignMismatchPenalty int
	swalignGapOpen         float32
	swalignGapExtend       float32
	swalignGapOpenIns      float32
	swalignGapExtendIns    float32
	swalignGapOpenDel      float32
	swalignGapExtendDel    float32
	swalignUseClipping     bool
	swalignClipOpen        float32
	swalignClipExtend      float32
	swalignHPOpenScale     float32
	swalignHPExtendScale   float32
	swalignHPOpenCap       float32
	swalignHPExtendCap     float32
	swalignGlobal          bool
	swalignVerbose         bool
)

func init() {
	swalignCmd.Flags().IntVar(&swalignMatchScore, "match", 1, "Match score")
	swalignCmd.Flags().IntVar(&swalignMismatchPenalty, "mismatch", 2, "Mismatch penalty")

	swalignCmd.Flags().Float32Var(&swalignGapOpenIns, "gap-open-ins", 6, "Insertion gap open penalty")
	swalignCmd.Flags().Float32Var(&swalignGapExtendIns, "gap-extend-ins", 1, "Insertion gap extension penalty")

	swalignCmd.Flags().Float32Var(&swalignGapOpenDel, "gap-open-del", 6, "Deletion gap open penalty")
	swalignCmd.Flags().Float32Var(&swalignGapExtendDel, "gap-extend-del", 1, "Deletion gap extension penalty")

	swalignCmd.Flags().Float32Var(&swalignGapOpen, "gap-open", 6, "Gap open penalty for both insertions and deletions (--gap-open-ins/-del override)")
	swalignCmd.Flags().Float32Var(&swalignGapExtend, "gap-extend", 1, "Gap extension penalty for both insertions and deletions (--gap-extend-ins/-del override)")

	swalignCmd.Flags().BoolVar(&swalignUseClipping, "clip", false, "Enable clipping penalties")
	swalignCmd.Flags().Float32Var(&swalignClipOpen, "clip-open", 5, "Clipping gap open penalty (used only when --clip is set)")
	swalignCmd.Flags().Float32Var(&swalignClipExtend, "clip-extend", 1, "Clipping gap extension penalty (used only when --clip is set)")

	swalignCmd.Flags().Float32Var(&swalignHPOpenScale, "hp-open-scale", 0, "Homopolymer gap-open discount scale")
	swalignCmd.Flags().Float32Var(&swalignHPOpenCap, "hp-open-cap", 0, "Homopolymer gap-open discount cap")
	swalignCmd.Flags().Float32Var(&swalignHPExtendScale, "hp-extend-scale", 0, "Homopolymer gap-extension discount scale")
	swalignCmd.Flags().Float32Var(&swalignHPExtendCap, "hp-extend-cap", 0, "Homopolymer gap-extension discount cap")
	swalignCmd.Flags().BoolVarP(&swalignVerbose, "verbose", "v", false, "Enable verbose aligner debug output")
	swalignCmd.Flags().BoolVar(&swalignGlobal, "global", false, "Enable global alignment (default local)")

}
