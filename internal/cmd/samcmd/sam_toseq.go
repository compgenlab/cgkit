package samcmd

import (
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/compgen-io/cgltk/htsio"
	"github.com/compgen-io/cgltk/seqio"
	"github.com/compgen-io/cgltk/support/sequtils"
	"github.com/compgen-io/cgltk/support/stringutils"
	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// tagFilter represents a single tag-based filter condition.
type tagFilter struct {
	tag string
	op  string // "eq", "not-eq", "contains", "not-contains", "lt", "gt", "lte", "gte"
	val string
}

// matches returns true if the SAM record passes this tag filter.
func (f *tagFilter) matches(rec *htsio.SamRecord) bool {
	t, ok := rec.Tags[f.tag]
	if !ok {
		return false
	}

	switch f.op {
	case "eq":
		return t.Value == f.val
	case "not-eq":
		return t.Value != f.val
	case "contains":
		return strings.Contains(t.Value, f.val)
	case "not-contains":
		return !strings.Contains(t.Value, f.val)
	case "lt", "gt", "lte", "gte":
		return f.numericCompare(t)
	}
	return false
}

func (f *tagFilter) numericCompare(t htsio.SamTag) bool {
	switch t.Type {
	case 'i':
		tv, ok := t.IntValue()
		if !ok {
			return false
		}
		fv, err := strconv.Atoi(f.val)
		if err != nil {
			return false
		}
		switch f.op {
		case "lt":
			return tv < fv
		case "gt":
			return tv > fv
		case "lte":
			return tv <= fv
		case "gte":
			return tv >= fv
		}
	case 'f':
		tv, ok := t.FloatValue()
		if !ok {
			return false
		}
		fv, err := strconv.ParseFloat(f.val, 64)
		if err != nil {
			return false
		}
		switch f.op {
		case "lt":
			return tv < fv
		case "gt":
			return tv > fv
		case "lte":
			return tv <= fv
		case "gte":
			return tv >= fv
		}
	}
	return false
}

// valStringArray wraps a pflag.Value to override Type() to "val".
type valStringArray struct {
	inner pflag.Value
}

func (v *valStringArray) String() string   { return v.inner.String() }
func (v *valStringArray) Set(s string) error { return v.inner.Set(s) }
func (v *valStringArray) Type() string     { return "val" }

// parseTagFilter parses "TAG:VALUE" into a tagFilter with the given op.
func parseTagFilter(s string, op string) (*tagFilter, error) {
	idx := strings.Index(s, ":")
	if idx < 1 {
		return nil, fmt.Errorf("invalid tag filter %q: expected TAG:VALUE", s)
	}
	return &tagFilter{
		tag: s[:idx],
		op:  op,
		val: s[idx+1:],
	}, nil
}

var samToSeqCmd = &cobra.Command{
	GroupID: "samcmd",
	Use:     "sam-toseq <input.bam> [output]",
	Short:   "Convert SAM/BAM/CRAM reads to FASTA or FASTQ",
	Long:    "Write SAM/BAM/CRAM reads to FASTA or FASTQ. Output file is optional; defaults to stdout.",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}

		if !samToSeqFasta && !samToSeqFastq {
			return fmt.Errorf("at least one of --fasta or --fastq must be specified")
		}

		// Parse tag filters.
		var filters []*tagFilter
		for _, spec := range samToSeqTagEq {
			f, err := parseTagFilter(spec, "eq")
			if err != nil {
				return err
			}
			filters = append(filters, f)
		}
		for _, spec := range samToSeqTagNotEq {
			f, err := parseTagFilter(spec, "not-eq")
			if err != nil {
				return err
			}
			filters = append(filters, f)
		}
		for _, spec := range samToSeqTagContains {
			f, err := parseTagFilter(spec, "contains")
			if err != nil {
				return err
			}
			filters = append(filters, f)
		}
		for _, spec := range samToSeqTagNotContains {
			f, err := parseTagFilter(spec, "not-contains")
			if err != nil {
				return err
			}
			filters = append(filters, f)
		}
		for _, spec := range samToSeqTagLt {
			f, err := parseTagFilter(spec, "lt")
			if err != nil {
				return err
			}
			filters = append(filters, f)
		}
		for _, spec := range samToSeqTagGt {
			f, err := parseTagFilter(spec, "gt")
			if err != nil {
				return err
			}
			filters = append(filters, f)
		}
		for _, spec := range samToSeqTagLte {
			f, err := parseTagFilter(spec, "lte")
			if err != nil {
				return err
			}
			filters = append(filters, f)
		}
		for _, spec := range samToSeqTagGte {
			f, err := parseTagFilter(spec, "gte")
			if err != nil {
				return err
			}
			filters = append(filters, f)
		}

		// Build SamReader options.
		opts := htsio.NewSamReaderOpts()
		if samToSeqFlagRequired != 0 {
			opts.FlagRequired(samToSeqFlagRequired)
		}
		if samToSeqFlagFilter != 0 {
			opts.FlagFilter(samToSeqFlagFilter)
		}
		if samToSeqMinMapQ > 0 {
			opts.MinMapQ(samToSeqMinMapQ)
		}
		if samToSeqRegion != "" {
			opts.Region(samToSeqRegion)
		} else if samToSeqRef != "" {
			opts.Region(samToSeqRef)
		}

		inputFile := args[0]
		reader, err := htsio.NewSamReader(inputFile, opts)
		if err != nil {
			return err
		}
		defer reader.Close()

		// Open output writers.
		var fastaWriter *seqio.FastaWriter
		var fastqWriter *seqio.FastqWriter

		if len(args) > 1 {
			outputFile := args[1]
			if samToSeqFasta {
				fastaWriter, err = seqio.OpenFastaWriter(outputFile)
				if err != nil {
					return fmt.Errorf("open FASTA output: %w", err)
				}
				defer fastaWriter.Close()
			} else {
				fastqWriter, err = seqio.OpenFastqWriter(outputFile)
				if err != nil {
					return fmt.Errorf("open FASTQ output: %w", err)
				}
				defer fastqWriter.Close()
			}
		} else {
			if samToSeqFasta {
				fastaWriter = seqio.NewFastaWriter(os.Stdout)
				defer fastaWriter.Close()
			} else {
				fastqWriter = seqio.NewFastqWriter(os.Stdout)
				defer fastqWriter.Close()
			}
		}

		for {
			rec, err := reader.Next()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}

			// Apply tag filters (ANDed).
			pass := true
			for _, f := range filters {
				if !f.matches(rec) {
					pass = false
					break
				}
			}
			if !pass {
				continue
			}

			seq := rec.Seq
			qual := rec.Qual

			// Reverse complement if on reverse strand.
			if rec.IsReverse() {
				seq = sequtils.ReverseCompliment(seq)
				qual = stringutils.ReverseString(qual)
			}

			if fastaWriter != nil {
				if err := fastaWriter.WriteRecord(rec.ReadName, "", seq); err != nil {
					return fmt.Errorf("write FASTA: %w", err)
				}
			}
			if fastqWriter != nil {
				if err := fastqWriter.WriteRecord(rec.ReadName, "", seq, qual); err != nil {
					return fmt.Errorf("write FASTQ: %w", err)
				}
			}
		}

		return nil
	},
}

var (
	samToSeqFasta        bool
	samToSeqFastq        bool
	samToSeqFlagRequired int
	samToSeqFlagFilter   int
	samToSeqMinMapQ      int
	samToSeqRegion       string
	samToSeqRef          string

	samToSeqTagEq          []string
	samToSeqTagNotEq       []string
	samToSeqTagContains    []string
	samToSeqTagNotContains []string
	samToSeqTagLt          []string
	samToSeqTagGt          []string
	samToSeqTagLte         []string
	samToSeqTagGte         []string
)

func init() {
	samToSeqCmd.Flags().BoolVar(&samToSeqFasta, "fasta", false, "Output in FASTA format")
	samToSeqCmd.Flags().BoolVar(&samToSeqFastq, "fastq", false, "Output in FASTQ format")

	samToSeqCmd.Flags().IntVar(&samToSeqFlagRequired, "flag-required", 0, "Require all of these flag bits to be set")
	samToSeqCmd.Flags().IntVar(&samToSeqFlagFilter, "flag-filter", 0, "Exclude reads with any of these flag bits set")
	samToSeqCmd.Flags().IntVar(&samToSeqMinMapQ, "min-mapq", 0, "Minimum mapping quality")
	samToSeqCmd.Flags().StringVar(&samToSeqRegion, "region", "", "Genomic region (chrom:start-end)")
	samToSeqCmd.Flags().StringVar(&samToSeqRef, "ref", "", "Filter by reference name")

	samToSeqCmd.Flags().StringArrayVar(&samToSeqTagEq, "tag-eq", nil, "Filter: tag equals value (TAG:VALUE)")
	samToSeqCmd.Flags().StringArrayVar(&samToSeqTagNotEq, "tag-not-eq", nil, "Filter: tag does not equal value (TAG:VALUE)")
	samToSeqCmd.Flags().StringArrayVar(&samToSeqTagContains, "tag-contains", nil, "Filter: tag contains substring (TAG:VALUE)")
	samToSeqCmd.Flags().StringArrayVar(&samToSeqTagNotContains, "tag-not-contains", nil, "Filter: tag does not contain substring (TAG:VALUE)")
	samToSeqCmd.Flags().StringArrayVar(&samToSeqTagLt, "tag-lt", nil, "Filter: tag less than value (TAG:VALUE)")
	samToSeqCmd.Flags().StringArrayVar(&samToSeqTagGt, "tag-gt", nil, "Filter: tag greater than value (TAG:VALUE)")
	samToSeqCmd.Flags().StringArrayVar(&samToSeqTagLte, "tag-lte", nil, "Filter: tag less than or equal to value (TAG:VALUE)")
	samToSeqCmd.Flags().StringArrayVar(&samToSeqTagGte, "tag-gte", nil, "Filter: tag greater than or equal to value (TAG:VALUE)")

	// Override the type placeholder from "stringArray" to "val" in help output.
	for _, name := range []string{"tag-eq", "tag-not-eq", "tag-contains", "tag-not-contains", "tag-lt", "tag-gt", "tag-lte", "tag-gte"} {
		samToSeqCmd.Flags().Lookup(name).Value = &valStringArray{inner: samToSeqCmd.Flags().Lookup(name).Value}
	}
}
