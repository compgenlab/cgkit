package bedcmd

import (
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/compgenlab/hts/bed"
	"github.com/spf13/cobra"
)

var (
	bedSetOutput       string
	bedSetInter        bool
	bedSetUnion        bool
	bedSetSub          bool
	bedSetExclusive    bool
	bedSetXor          bool
	bedSetIgnoreStrand bool
	bedSetBed3         bool
	bedSetSum          bool
	bedSetCount        bool
	bedSetDelim        string
	bedSetLabelA       string
	bedSetLabelB       string
	bedSetFlanking     int
	bedSetTabix        bool
)

var bedSetCmd = &cobra.Command{
	GroupID:     "bedcmd",
	Annotations: map[string]string{"since": "v0.3.1"},
	Use:         "bed-set <mode> <A.bed> <B.bed>",
	Short:       "Set algebra (intersect/union/subtract/exclusive) on two BED files",
	Long: `Coordinate set algebra over two BED files. Exactly one mode is required:

  --inter       bases covered by both A and B
  --union       bases covered by A or B (merged)
  --sub         bases in A but not B (A - B)
  --exclusive   bases in exactly one of A or B (alias: --xor)

Inputs should be coordinate-sorted; output is emitted sorted. When both inputs
are BED6 the operation is strand-aware (regions overlap only on the same strand)
and output is BED6; --ignore-strand (or --bed3) collapses strand and emits BED3.

In BED6 output the name column is the merged contributing names (--delim,
default "|") or provenance labels (--a/--b), and the score column is 0 unless
--sum or --count is given. -o NAME.gz writes sorted BGZF; --tbi adds a tabix index.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) < 2 {
			cmd.Help()
			return nil
		}

		// Resolve the single mode.
		op := bed.OpInter
		modes := 0
		if bedSetInter {
			modes++
			op = bed.OpInter
		}
		if bedSetUnion {
			modes++
			op = bed.OpUnion
		}
		if bedSetSub {
			modes++
			op = bed.OpSub
		}
		if bedSetExclusive || bedSetXor {
			modes++
			op = bed.OpXor
		}
		if modes != 1 {
			return fmt.Errorf("exactly one of --inter, --union, --sub, --exclusive is required")
		}
		if bedSetSum && bedSetCount {
			return fmt.Errorf("--sum and --count are mutually exclusive")
		}
		if bedSetFlanking != 0 && op != bed.OpUnion {
			return fmt.Errorf("--flanking is only valid with --union")
		}
		if args[0] == "-" && args[1] == "-" {
			return fmt.Errorf("only one input may be read from stdin")
		}

		aRecs, aCols, err := readBedAll(cmd, args[0])
		if err != nil {
			return err
		}
		bRecs, bCols, err := readBedAll(cmd, args[1])
		if err != nil {
			return err
		}

		aStranded := aCols >= 6
		bStranded := bCols >= 6
		strandAware := !bedSetIgnoreStrand && !bedSetBed3 && aStranded && bStranded
		if !bedSetIgnoreStrand && !bedSetBed3 && aStranded != bStranded {
			fmt.Fprintln(os.Stderr, "warning: inputs have mixed column widths; treating as strand-agnostic (BED3). Use --ignore-strand to silence.")
		}

		writer, err := openBedSetWriter(cmd, strandAware)
		if err != nil {
			return err
		}

		opts := bed.SetOpts{Op: op, IgnoreStrand: !strandAware, FlankX: bedSetFlanking}
		for seg := range bed.SetOperationRecords(aRecs, bRecs, opts) {
			rec := &bed.BedRecord{
				Ref:     seg.Ref,
				Start:   seg.Start,
				End:     seg.End,
				Strand:  seg.Strand,
				HasName: true,
				Name:    bedSetSegName(seg),
				Score:   bedSetSegScore(seg),
			}
			if err := writer.WriteRecord(rec); err != nil {
				return err
			}
		}
		return writer.Close()
	},
}

// readBedAll reads every record from filename ("-" = stdin) and reports the
// maximum column count seen (to detect BED6 inputs).
func readBedAll(cmd *cobra.Command, filename string) ([]*bed.BedRecord, int, error) {
	var r *bed.BedReader
	var err error
	if filename == "-" {
		r, err = bed.NewBedReader(cmd.InOrStdin())
	} else {
		r, err = bed.NewBedFile(filename)
	}
	if err != nil {
		return nil, 0, err
	}
	defer r.Close()

	var recs []*bed.BedRecord
	maxCols := 0
	for {
		rec, err := r.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, 0, err
		}
		if rec.Columns > maxCols {
			maxCols = rec.Columns
		}
		recs = append(recs, rec)
	}
	return recs, maxCols, nil
}

// openBedSetWriter selects the output sink and column layout.
func openBedSetWriter(cmd *cobra.Command, strandAware bool) (*bed.BedWriter, error) {
	opts := bed.NewBedWriterOpts()
	if strandAware {
		opts.Columns(bed.Columns6).RawStrand(true)
	} else {
		opts.Columns(bed.Columns3)
	}
	if bedSetCount {
		opts.ForceScoreInt(true)
	}

	if bedSetTabix {
		if bedSetOutput == "" || bedSetOutput == "-" {
			return nil, fmt.Errorf("--tbi requires an output file (-o)")
		}
		return bed.OpenBedWriter(bedSetOutput, opts.Tabix(true))
	}
	if bedSetOutput == "" || bedSetOutput == "-" {
		return bed.NewBedWriter(cmd.OutOrStdout(), opts), nil
	}
	if strings.HasSuffix(bedSetOutput, ".gz") {
		return bed.OpenBedWriter(bedSetOutput, opts.Bgzip(true))
	}
	return bed.OpenBedWriter(bedSetOutput, opts)
}

// bedSetSegName builds the col-4 name: provenance labels when --a/--b are set,
// otherwise the deduped merged contributing names. "." when empty.
func bedSetSegName(seg *bed.OutSegment) string {
	if bedSetLabelA != "" || bedSetLabelB != "" {
		hasA, hasB := false, false
		for _, c := range seg.Contribs {
			if c.Source == 0 {
				hasA = true
			} else {
				hasB = true
			}
		}
		var parts []string
		if hasA && bedSetLabelA != "" {
			parts = append(parts, bedSetLabelA)
		}
		if hasB && bedSetLabelB != "" {
			parts = append(parts, bedSetLabelB)
		}
		if len(parts) == 0 {
			return "."
		}
		return strings.Join(parts, bedSetDelim)
	}

	var parts []string
	seen := map[string]bool{}
	for _, c := range seg.Contribs {
		if c.Name != "" && !seen[c.Name] {
			seen[c.Name] = true
			parts = append(parts, c.Name)
		}
	}
	if len(parts) == 0 {
		return "."
	}
	return strings.Join(parts, bedSetDelim)
}

// bedSetSegScore builds the col-5 score: sum of contributing scores, count of
// contributing regions, or 0.
func bedSetSegScore(seg *bed.OutSegment) float64 {
	switch {
	case bedSetSum:
		var sum float64
		for _, c := range seg.Contribs {
			sum += c.Score
		}
		return sum
	case bedSetCount:
		return float64(len(seg.Contribs))
	default:
		return 0
	}
}

func init() {
	bedSetCmd.Flags().BoolVar(&bedSetInter, "inter", false, "Intersection: bases covered by both A and B")
	bedSetCmd.Flags().BoolVar(&bedSetUnion, "union", false, "Union: bases covered by A or B (merged)")
	bedSetCmd.Flags().BoolVar(&bedSetSub, "sub", false, "Subtraction: bases in A but not B (A - B)")
	bedSetCmd.Flags().BoolVar(&bedSetExclusive, "exclusive", false, "Exclusive: bases in exactly one of A or B")
	bedSetCmd.Flags().BoolVar(&bedSetXor, "xor", false, "Alias for --exclusive")
	bedSetCmd.Flags().BoolVar(&bedSetIgnoreStrand, "ignore-strand", false, "Ignore strand; emit BED3")
	bedSetCmd.Flags().BoolVar(&bedSetBed3, "bed3", false, "Force bare BED3 output")
	bedSetCmd.Flags().BoolVar(&bedSetSum, "sum", false, "BED6 score = sum of contributing region scores")
	bedSetCmd.Flags().BoolVar(&bedSetCount, "count", false, "BED6 score = number of contributing regions")
	bedSetCmd.Flags().StringVar(&bedSetDelim, "delim", "|", "Delimiter for merged names / provenance labels")
	bedSetCmd.Flags().StringVar(&bedSetLabelA, "a", "", "Provenance label for regions from A (sets name by source)")
	bedSetCmd.Flags().StringVar(&bedSetLabelB, "b", "", "Provenance label for regions from B (sets name by source)")
	bedSetCmd.Flags().IntVar(&bedSetFlanking, "flanking", 0, "Union only: merge regions within N bases of each other")
	bedSetCmd.Flags().BoolVar(&bedSetTabix, "tbi", false, "Write a tabix (.tbi) index (implies sorted BGZF; requires -o file)")
	bedSetCmd.Flags().StringVarP(&bedSetOutput, "output", "o", "-", "Output filename (.gz = sorted BGZF; - for stdout)")
}
