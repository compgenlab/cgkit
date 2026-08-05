package bedcmd

import (
	"fmt"
	"io"
	"sort"
	"strconv"

	"github.com/compgenlab/cgkit/internal/cmdio"
	"github.com/spf13/cobra"
)

var bedStatsOutput string

var bedStatsCmd = &cobra.Command{
	GroupID:     "bedcmd",
	Annotations: map[string]string{"since": "v0.3.1"},
	Use:         "bed-stats <input.bed>",
	Short:       "Summary statistics for a BED file",
	Args:        cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {

		reader, err := openBedInput(cmd, args[0])
		if err != nil {
			return err
		}
		defer reader.Close()

		var total, totalSize int
		sizes := newTallyCounts()
		refOrder := []string{}
		refSeen := map[string]bool{}
		refCount := map[string]int{}

		for {
			rec, err := reader.NextRecord()
			if err != nil {
				if err == io.EOF {
					break
				}
				return err
			}

			length := rec.Length()
			total++
			totalSize += length
			sizes.incr(length)

			if !refSeen[rec.Ref] {
				refSeen[rec.Ref] = true
				refOrder = append(refOrder, rec.Ref)
			}
			refCount[rec.Ref]++
		}

		out := cmd.OutOrStdout()
		dst, err := cmdio.CreateOutput(cmd, "-o/--output", bedStatsOutput)
		if err != nil {
			return err
		}
		// Deferred for the error paths, closed explicitly below so the close
		// error reaches the exit status. Output.Close is safe twice.
		defer dst.Close()
		out = dst.W

		fmt.Fprintf(out, "Total number of regions:\t%d\n", total)
		fmt.Fprintf(out, "Total number of bases:\t%d\n", totalSize)
		fmt.Fprintln(out, "")
		fmt.Fprintf(out, "Mean region size:\t%s\n", formatMean(sizes.mean()))
		fmt.Fprintf(out, "Median region size:\t%d\n", sizes.median())
		fmt.Fprintf(out, "Max region size:\t%d\n", sizes.max)
		fmt.Fprintf(out, "Min region size:\t%d\n", sizes.min)
		fmt.Fprintln(out, "")
		for _, ref := range refOrder {
			fmt.Fprintf(out, "%s\t%d\n", ref, refCount[ref])
		}
		return dst.Close()
	},
}

// tallyCounts mirrors the subset of io.compgen.common.TallyCounts used by
// bed-stats: negative keys are ignored, min/max start at -1, and the mean and
// median are computed over the tallied (non-negative) values.
type tallyCounts struct {
	counts     map[int]int
	min        int
	max        int
	sum        int
	totalCount int
}

func newTallyCounts() *tallyCounts {
	return &tallyCounts{counts: map[int]int{}, min: -1, max: -1}
}

func (t *tallyCounts) incr(k int) {
	if k < 0 {
		return
	}
	t.totalCount++
	t.sum += k
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

// mean returns the arithmetic mean of the tallied values, or NaN when empty
// (matching Java's (double)0/0).
func (t *tallyCounts) mean() float64 {
	return float64(t.sum) / float64(t.totalCount)
}

// median returns the 0.5 quantile using the same accumulation rule as
// TallyCounts.getQuantile: walk values in ascending order until the cumulative
// count exceeds 0.5*totalCount, else return max.
func (t *tallyCounts) median() int {
	thres := 0.5 * float64(t.totalCount)
	keys := make([]int, 0, len(t.counts))
	for k := range t.counts {
		keys = append(keys, k)
	}
	sort.Ints(keys)

	count := 0
	for _, k := range keys {
		if !(float64(count) < thres) {
			break
		}
		count += t.counts[k]
		if float64(count) > thres {
			return k
		}
	}
	return t.max
}

// formatMean renders the mean as a plain decimal (shortest round-trip). Unlike
// Java's Double.toString it does not emit scientific notation or a trailing
// ".0" for integral means; see the project notes on score formatting.
func formatMean(f float64) string {
	return strconv.FormatFloat(f, 'g', -1, 64)
}

func init() {
	bedStatsCmd.Flags().StringVarP(&bedStatsOutput, "output", "o", "-", "Output filename (- for stdout)")
}
