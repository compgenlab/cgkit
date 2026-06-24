package vcfcmd

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strings"

	"github.com/spf13/cobra"
)

var (
	vcfStatsPassing     bool
	vcfStatsFilterCombo bool
	vcfStatsInfoTally   []string
	vcfStatsInfoPresent []string
	vcfStatsRegion      string
)

var vcfStatsCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-stats <input.vcf>",
	Short:       "Summary statistics about a VCF file",
	Long: `Summary statistics about a VCF file: variant counts, SNV/indel/reference-only
breakdown, Ts/Tv counts and ratio, and filter tallies.

--info-tally counts the distinct values of an INFO field; --info-present counts
how often a field is present or absent. Both may be repeated or
comma-separated.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		src, err := openRecordSource(cmd, args[0], vcfStatsRegion)
		if err != nil {
			return err
		}
		defer src.close()

		infoFields := splitAll(vcfStatsInfoTally)
		presentFields := splitAll(vcfStatsInfoPresent)

		var count, passing, filtered, refonly, tsCount, tvCount, indel int64
		filterCounts := map[string]int64{}
		fullFilterCounts := map[string]int64{}
		infoTally := make([]map[string]int64, len(infoFields))
		infoMissing := make([]int64, len(infoFields))
		for i := range infoTally {
			infoTally[i] = map[string]int64{}
		}
		present := make([]int64, len(presentFields))
		absent := make([]int64, len(presentFields))

		for {
			rec, err := src.next()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			count++

			if rec.IsFiltered() {
				filtered++
				for _, f := range rec.Filters() {
					filterCounts[f]++
				}
				if vcfStatsFilterCombo {
					filters := append([]string(nil), rec.Filters()...)
					sort.Strings(filters)
					fullFilterCounts[strings.Join(filters, ",")]++
				}
				if vcfStatsPassing {
					continue
				}
			} else {
				passing++
			}

			if len(rec.Alt()) == 0 {
				refonly++
			} else if rec.IsIndel() {
				indel++
			}

			switch rec.CalcTsTv() {
			case -1:
				tsCount++
			case 1:
				tvCount++
			}

			for i, field := range infoFields {
				if v, ok := rec.Info().Get(field); ok {
					infoTally[i][v.String()]++
				} else {
					infoMissing[i]++
				}
			}
			for i, field := range presentFields {
				if rec.Info().Contains(field) {
					present[i]++
				} else {
					absent[i]++
				}
			}
		}

		out := bufio.NewWriter(cmd.OutOrStdout())
		fmt.Fprintf(out, "Total variants:\t%d\n", count)
		fmt.Fprintf(out, "Filtered variants:\t%d\n", filtered)
		fmt.Fprintf(out, "Passing variants:\t%d\n", passing)
		fmt.Fprintln(out)
		fmt.Fprintf(out, "SNV:\t%d\n", count-indel-refonly)
		fmt.Fprintf(out, "Indels:\t%d\n", indel)
		fmt.Fprintf(out, "Reference-only:\t%d\n", refonly)
		fmt.Fprintln(out)
		fmt.Fprintf(out, "Transitions:\t%d\n", tsCount)
		fmt.Fprintf(out, "Transversions:\t%d\n", tvCount)
		ratio := ""
		if tvCount > 0 {
			ratio = javaDouble(float64(tsCount) / float64(tvCount))
		}
		fmt.Fprintf(out, "Ts/Tv ratio:\t%s\n", ratio)
		fmt.Fprintln(out)

		if vcfStatsFilterCombo {
			fmt.Fprintln(out)
			fmt.Fprintln(out, "[Filter combinations]")
			for _, k := range sortedKeys(fullFilterCounts) {
				fmt.Fprintf(out, "%s: %d\n", k, fullFilterCounts[k])
			}
		} else {
			fmt.Fprintln(out, "[Filters]")
			for _, k := range sortedKeys(filterCounts) {
				fmt.Fprintf(out, "%s: %d\n", k, filterCounts[k])
			}
		}

		for i, field := range infoFields {
			fmt.Fprintln(out)
			fmt.Fprintf(out, "[%s]\n", field)
			for _, k := range sortedKeys(infoTally[i]) {
				fmt.Fprintf(out, "%s\t%d\n", k, infoTally[i][k])
			}
			fmt.Fprintf(out, "*missing*\t%d\n", infoMissing[i])
		}
		for i, field := range presentFields {
			fmt.Fprintln(out)
			fmt.Fprintf(out, "[%s]\n", field)
			fmt.Fprintf(out, "Present\t%d\n", present[i])
			fmt.Fprintf(out, "Absent\t%d\n", absent[i])
		}

		return out.Flush()
	},
}

// splitAll splits each comma-separated value and flattens the result.
func splitAll(vals []string) []string {
	var out []string
	for _, v := range vals {
		out = append(out, strings.Split(v, ",")...)
	}
	return out
}

func sortedKeys(m map[string]int64) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

func init() {
	vcfStatsCmd.Flags().BoolVar(&vcfStatsPassing, "passing", false, "Only process passing variants")
	vcfStatsCmd.Flags().BoolVar(&vcfStatsFilterCombo, "filter-combo", false, "Show full filter combinations")
	vcfStatsCmd.Flags().StringArrayVar(&vcfStatsInfoTally, "info-tally", nil, "Tally the values of these INFO fields (repeatable, comma-separated)")
	vcfStatsCmd.Flags().StringArrayVar(&vcfStatsInfoPresent, "info-present", nil, "Tally presence/absence of these INFO fields (repeatable, comma-separated)")
	vcfStatsCmd.Flags().StringVar(&vcfStatsRegion, "region", "", "Only variants in this 1-based region (chrom:start-end, or chrom); requires a tabix-indexed file")
}
