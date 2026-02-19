package cmd

import (
	"fmt"

	"github.com/compgen-io/cgltk/align"
	"github.com/compgen-io/cgltk/sequtils"
	"github.com/spf13/cobra"
)

// fastagcCmd implements the initial counting entrypoint.
var swalignCmd = &cobra.Command{
	Use:   "pairwise query target",
	Short: "Align the two given sequences",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) != 2 {
			cmd.Help()
			return nil
		}

		if swalignUseClipping && swalignGlobal {
			return fmt.Errorf("cannot use clipping penalties with global alignment")
		}

		opts := align.DnaAlignmentDefaults().
			ScoringMatrix(align.MatchMismatchScoring(swalignMatchScore, swalignMismatchPenalty)).
			GapPenaltyIns(swalignGapOpenIns, swalignGapExtendIns).
			GapPenaltyDel(swalignGapOpenDel, swalignGapExtendDel).
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
		aln1 := sw.Align(args[0], args[1])
		aln2 := sw.Align(sequtils.ReverseCompliment(args[0]), args[1])
		if aln1.Score >= aln2.Score {
			fmt.Println(aln1.String())
		} else {
			aln2.QueryRevComp = true
			fmt.Println(aln2.String())
		}
		return nil
	},
}

var (
	swalignMatchScore      int
	swalignMismatchPenalty int
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

	swalignCmd.Flags().Float32Var(&swalignGapOpenIns, "gap-open-ins", 5, "Insertion gap open penalty")
	swalignCmd.Flags().Float32Var(&swalignGapExtendIns, "gap-extend-ins", 2, "Insertion gap extension penalty")

	swalignCmd.Flags().Float32Var(&swalignGapOpenDel, "gap-open-del", 5, "Deletion gap open penalty")
	swalignCmd.Flags().Float32Var(&swalignGapExtendDel, "gap-extend-del", 2, "Deletion gap extension penalty")

	swalignCmd.Flags().BoolVar(&swalignUseClipping, "clip", false, "Enable clipping penalties")
	swalignCmd.Flags().Float32Var(&swalignClipOpen, "clip-open", 5, "Clipping gap open penalty (used only when --clip is set)")
	swalignCmd.Flags().Float32Var(&swalignClipExtend, "clip-extend", 1, "Clipping gap extension penalty (used only when --clip is set)")

	swalignCmd.Flags().Float32Var(&swalignHPOpenScale, "hp-open-scale", 0, "Homopolymer gap-open discount scale")
	swalignCmd.Flags().Float32Var(&swalignHPOpenCap, "hp-open-cap", 0, "Homopolymer gap-open discount cap")
	swalignCmd.Flags().Float32Var(&swalignHPExtendScale, "hp-extend-scale", 0, "Homopolymer gap-extension discount scale")
	swalignCmd.Flags().Float32Var(&swalignHPExtendCap, "hp-extend-cap", 0, "Homopolymer gap-extension discount cap")
	swalignCmd.Flags().BoolVarP(&swalignVerbose, "verbose", "v", false, "Enable verbose aligner debug output")
	swalignCmd.Flags().BoolVarP(&swalignGlobal, "global", "g", false, "Enable global alignment")

	rootCmd.AddCommand(swalignCmd)
}
