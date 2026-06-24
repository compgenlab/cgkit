package vcfcmd

import (
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/compgenlab/hts/vcf"
	"github.com/compgenlab/hts/vcf/filter"
	"github.com/spf13/cobra"
)

var (
	vcfFilterOutput  string
	vcfFilterPassing bool
	vcfFilterFailing bool
	vcfFilterStats   string
	vcfFilterChain   []chainArg
)

var vcfFilterCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-filter <input.vcf>",
	Short:       "Filter a VCF file by stamping FILTER codes",
	Long: `Filter a VCF file. Each flag adds a filter; every variant runs through the
filters in command-line order, and a variant that fails a filter gets that
filter's code stamped into its FILTER column. A variant with no codes is PASS.

  --chrom-pass CSV   flag variants not on these chromosomes
  --chrom-fail CSV   flag variants on these chromosomes
  --snv / --indel    flag SNVs / indels
  --max-ins N        flag insertions longer than N
  --max-del N        flag deletions longer than N
  --qual FLOAT       flag variants with QUAL below FLOAT
  --het / --hom      flag heterozygous / homozygous variants (requires GT)
  --flag-present KEY  flag variants whose INFO contains the flag KEY
  --flag-missing KEY  flag variants whose INFO is missing the flag KEY
  --value-missing KEY[:SAMPLE]    flag variants missing a value
  --eq / --neq KEY:VAL[:SAMPLE[:ALLELE]]          string (in)equality
  --contain / --not-contain KEY:VAL[:SAMPLE[:ALLELE]]   substring match
  --in / --not-in KEY:CSV[:SAMPLE[:ALLELE]]       membership in a CSV list
  --lt / --lte / --gt / --gte KEY:VAL[:SAMPLE[:ALLELE]]  numeric comparison

For the value filters, SAMPLE selects a sample (or INFO for an INFO field; empty
means every sample), and ALLELE selects an allele (ref, alt1, an index, or sum).

By default every variant is written (with FILTER updated). --passing writes only
variants that pass all filters; --failing writes only variants that fail one.
--stats FILE writes per-filter-combination counts (tab-separated).`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		chain, err := buildFilterChain()
		if err != nil {
			return err
		}

		reader, err := openVcfInput(cmd, args[0])
		if err != nil {
			return err
		}
		defer reader.Close()

		header, err := reader.Header()
		if err != nil {
			return err
		}
		if err := chain.SetupHeaders(header); err != nil {
			return err
		}
		stampVcfProvenance(header, "vcf-filter")

		var writer *vcf.VcfWriter
		var closeFile func() error
		if vcfFilterOutput == "" || vcfFilterOutput == "-" {
			writer = vcf.NewVcfWriter(cmd.OutOrStdout())
		} else {
			w, err := vcf.OpenVcfWriter(vcfFilterOutput)
			if err != nil {
				return err
			}
			writer = w
			closeFile = w.Close
		}
		if err := writer.WriteHeader(header); err != nil {
			return err
		}

		stats := newFilterStats()
		for {
			rec, err := reader.NextRecord()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if err := chain.Apply(rec); err != nil {
				return err
			}
			filtered := rec.IsFiltered()
			if filtered {
				stats.tally(rec.Filters())
			}
			if vcfFilterPassing && filtered {
				continue
			}
			if vcfFilterFailing && !filtered {
				continue
			}
			if err := writer.WriteRecord(rec); err != nil {
				return err
			}
		}
		if err := chain.Close(); err != nil {
			return err
		}
		if vcfFilterStats != "" {
			if err := stats.write(vcfFilterStats); err != nil {
				return err
			}
		}
		if closeFile != nil {
			return closeFile()
		}
		return writer.Close()
	},
}

// buildFilterChain assembles the filter chain from the flags in command-line
// order.
func buildFilterChain() (*filter.Chain, error) {
	c := filter.NewChain()
	for _, fl := range vcfFilterChain {
		switch fl.kind {
		case "snv":
			c.Add(filter.NewSNV())
		case "indel":
			c.Add(filter.NewIndel())
		case "het":
			c.Add(filter.NewHeterozygous())
		case "hom":
			c.Add(filter.NewHomozygous())
		case "chrom-pass":
			c.Add(filter.NewChromPass(strings.Split(fl.value, ",")))
		case "chrom-fail":
			c.Add(filter.NewChromFail(strings.Split(fl.value, ",")))
		case "max-ins":
			n, err := strconv.Atoi(fl.value)
			if err != nil {
				return nil, fmt.Errorf("--max-ins: %w", err)
			}
			c.Add(filter.NewMaxIns(n))
		case "max-del":
			n, err := strconv.Atoi(fl.value)
			if err != nil {
				return nil, fmt.Errorf("--max-del: %w", err)
			}
			c.Add(filter.NewMaxDel(n))
		case "qual":
			q, err := strconv.ParseFloat(fl.value, 64)
			if err != nil {
				return nil, fmt.Errorf("--qual: %w", err)
			}
			c.Add(filter.NewQual(q))
		case "flag-present":
			c.Add(filter.NewFlagPresent(fl.value))
		case "flag-missing":
			c.Add(filter.NewFlagAbsent(fl.value))
		case "value-missing":
			key, sample, err := splitMissingArg(fl.value)
			if err != nil {
				return nil, err
			}
			c.Add(filter.NewValueMissing(key, sample))
		case "eq", "neq", "contain", "not-contain", "in", "not-in", "lt", "lte", "gt", "gte":
			key, val, sample, allele, err := splitFilterArg(fl.kind, fl.value)
			if err != nil {
				return nil, err
			}
			f, err := newValueFilter(fl.kind, key, val, sample, allele)
			if err != nil {
				return nil, err
			}
			c.Add(f)
		}
	}
	return c, nil
}

// filterStats tallies, per distinct combination of FILTER codes, how many failing
// variants carried exactly that combination (ports ngsutilsj's --stats output).
type filterStats struct {
	order  []string
	counts map[string]int
}

func newFilterStats() *filterStats {
	return &filterStats{counts: map[string]int{}}
}

func (s *filterStats) tally(codes []string) {
	sorted := append([]string(nil), codes...)
	sort.Strings(sorted)
	key := strings.Join(sorted, ",")
	if _, ok := s.counts[key]; !ok {
		s.order = append(s.order, key)
	}
	s.counts[key]++
}

func (s *filterStats) write(filename string) error {
	out, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer out.Close()
	for _, key := range s.order {
		if _, err := fmt.Fprintf(out, "%s\t%d\n", key, s.counts[key]); err != nil {
			return err
		}
	}
	return nil
}

// splitFilterArg parses a KEY:VAL[:SAMPLE[:ALLELE]] filter argument.
func splitFilterArg(name, arg string) (key, val, sample, allele string, err error) {
	spl := strings.Split(arg, ":")
	switch len(spl) {
	case 2:
		return spl[0], spl[1], "", "", nil
	case 3:
		return spl[0], spl[1], spl[2], "", nil
	case 4:
		return spl[0], spl[1], spl[2], spl[3], nil
	default:
		return "", "", "", "", fmt.Errorf("--%s: malformed argument %q (want KEY:VAL[:SAMPLE[:ALLELE]])", name, arg)
	}
}

// splitMissingArg parses a KEY[:SAMPLE] --value-missing argument.
func splitMissingArg(arg string) (key, sample string, err error) {
	spl := strings.Split(arg, ":")
	switch len(spl) {
	case 1:
		return spl[0], "", nil
	case 2:
		return spl[0], spl[1], nil
	default:
		return "", "", fmt.Errorf("--value-missing: malformed argument %q (want KEY[:SAMPLE])", arg)
	}
}

// newValueFilter builds the comparison filter for the given flag kind.
func newValueFilter(kind, key, val, sample, allele string) (filter.Filter, error) {
	switch kind {
	case "eq":
		return filter.NewEquals(key, val, sample, allele), nil
	case "neq":
		return filter.NewNotEquals(key, val, sample, allele), nil
	case "contain":
		return filter.NewContains(key, val, sample, allele), nil
	case "not-contain":
		return filter.NewNotContains(key, val, sample, allele), nil
	case "in":
		return filter.NewInList(key, strings.Split(val, ","), sample, allele), nil
	case "not-in":
		return filter.NewNotInList(key, strings.Split(val, ","), sample, allele), nil
	case "lt", "lte", "gt", "gte":
		thres, err := strconv.ParseFloat(val, 64)
		if err != nil {
			return nil, fmt.Errorf("--%s: %q is not a number", kind, val)
		}
		switch kind {
		case "lt":
			return filter.NewLessThan(key, thres, sample, allele), nil
		case "lte":
			return filter.NewLessThanEqual(key, thres, sample, allele), nil
		case "gt":
			return filter.NewGreaterThan(key, thres, sample, allele), nil
		default:
			return filter.NewGreaterThanEqual(key, thres, sample, allele), nil
		}
	}
	return nil, fmt.Errorf("unknown filter kind %q", kind)
}

func init() {
	f := vcfFilterCmd.Flags()
	f.StringVarP(&vcfFilterOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
	f.BoolVar(&vcfFilterPassing, "passing", false, "Only output passing variants")
	f.BoolVar(&vcfFilterFailing, "failing", false, "Only output failing variants")
	f.StringVar(&vcfFilterStats, "stats", "", "Write per-filter-combination counts (tab-separated) to this file")
	registerChainVal(f, &vcfFilterChain, "chrom-pass", "Flag variants not on these chromosomes (CSV)")
	registerChainVal(f, &vcfFilterChain, "chrom-fail", "Flag variants on these chromosomes (CSV)")
	registerChainBool(f, &vcfFilterChain, "snv", "Flag SNVs")
	registerChainBool(f, &vcfFilterChain, "indel", "Flag indels")
	registerChainVal(f, &vcfFilterChain, "max-ins", "Flag insertions longer than this length")
	registerChainVal(f, &vcfFilterChain, "max-del", "Flag deletions longer than this length")
	registerChainVal(f, &vcfFilterChain, "qual", "Flag variants with QUAL below this value")
	registerChainBool(f, &vcfFilterChain, "het", "Flag heterozygous variants (requires GT)")
	registerChainBool(f, &vcfFilterChain, "hom", "Flag homozygous variants (requires GT)")
	registerChainVal(f, &vcfFilterChain, "flag-present", "Flag variants whose INFO contains this flag (KEY)")
	registerChainVal(f, &vcfFilterChain, "flag-missing", "Flag variants whose INFO is missing this flag (KEY)")
	registerChainVal(f, &vcfFilterChain, "value-missing", "Flag variants missing a value (KEY[:SAMPLE], SAMPLE=INFO for an INFO field)")
	registerChainVal(f, &vcfFilterChain, "eq", "Flag variants where KEY equals VAL (KEY:VAL[:SAMPLE[:ALLELE]])")
	registerChainVal(f, &vcfFilterChain, "neq", "Flag variants where KEY does not equal VAL (KEY:VAL[:SAMPLE[:ALLELE]])")
	registerChainVal(f, &vcfFilterChain, "contain", "Flag variants where KEY contains VAL (KEY:VAL[:SAMPLE[:ALLELE]])")
	registerChainVal(f, &vcfFilterChain, "not-contain", "Flag variants where KEY does not contain VAL (KEY:VAL[:SAMPLE[:ALLELE]])")
	registerChainVal(f, &vcfFilterChain, "in", "Flag variants where KEY is in the CSV list VAL (KEY:VAL[:SAMPLE[:ALLELE]])")
	registerChainVal(f, &vcfFilterChain, "not-in", "Flag variants where KEY is not in the CSV list VAL (KEY:VAL[:SAMPLE[:ALLELE]])")
	registerChainVal(f, &vcfFilterChain, "lt", "Flag variants where KEY < VAL (KEY:VAL[:SAMPLE[:ALLELE]])")
	registerChainVal(f, &vcfFilterChain, "lte", "Flag variants where KEY <= VAL (KEY:VAL[:SAMPLE[:ALLELE]])")
	registerChainVal(f, &vcfFilterChain, "gt", "Flag variants where KEY > VAL (KEY:VAL[:SAMPLE[:ALLELE]])")
	registerChainVal(f, &vcfFilterChain, "gte", "Flag variants where KEY >= VAL (KEY:VAL[:SAMPLE[:ALLELE]])")
}
