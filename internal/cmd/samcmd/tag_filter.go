package samcmd

import (
	"github.com/compgen-io/cgltk/htsio"
	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// valStringArray wraps a pflag.Value to override Type() to "val".
type valStringArray struct {
	inner pflag.Value
}

func (v *valStringArray) String() string     { return v.inner.String() }
func (v *valStringArray) Set(s string) error { return v.inner.Set(s) }
func (v *valStringArray) Type() string       { return "val" }

// tagFilterFlags holds the raw string slices from cobra flags.
type tagFilterFlags struct {
	Eq, NotEq, Contains, NotContains *[]string
	Lt, Gt, Lte, Gte                *[]string
}

// register adds tag filter flags to the command.
func (t *tagFilterFlags) register(cmd *cobra.Command) {
	t.Eq = new([]string)
	t.NotEq = new([]string)
	t.Contains = new([]string)
	t.NotContains = new([]string)
	t.Lt = new([]string)
	t.Gt = new([]string)
	t.Lte = new([]string)
	t.Gte = new([]string)

	cmd.Flags().StringArrayVar(t.Eq, "tag-eq", nil, "Keep reads where tag equals value (TAG:VALUE)")
	cmd.Flags().StringArrayVar(t.NotEq, "tag-not-eq", nil, "Keep reads where tag does not equal value (TAG:VALUE)")
	cmd.Flags().StringArrayVar(t.Contains, "tag-contains", nil, "Keep reads where tag contains substring (TAG:VALUE)")
	cmd.Flags().StringArrayVar(t.NotContains, "tag-not-contains", nil, "Keep reads where tag does not contain substring (TAG:VALUE)")
	cmd.Flags().StringArrayVar(t.Lt, "tag-lt", nil, "Keep reads where tag is less than value (TAG:VALUE)")
	cmd.Flags().StringArrayVar(t.Gt, "tag-gt", nil, "Keep reads where tag is greater than value (TAG:VALUE)")
	cmd.Flags().StringArrayVar(t.Lte, "tag-lte", nil, "Keep reads where tag is less than or equal to value (TAG:VALUE)")
	cmd.Flags().StringArrayVar(t.Gte, "tag-gte", nil, "Keep reads where tag is greater than or equal to value (TAG:VALUE)")

	for _, name := range []string{"tag-eq", "tag-not-eq", "tag-contains", "tag-not-contains", "tag-lt", "tag-gt", "tag-lte", "tag-gte"} {
		cmd.Flags().Lookup(name).Value = &valStringArray{inner: cmd.Flags().Lookup(name).Value}
	}
}

// parse converts the raw flag strings into htsio.TagFilter values.
func (t *tagFilterFlags) parse() ([]*htsio.TagFilter, error) {
	var filters []*htsio.TagFilter
	for _, pair := range []struct {
		specs *[]string
		op    htsio.TagFilterOp
	}{
		{t.Eq, htsio.TagEq}, {t.NotEq, htsio.TagNotEq},
		{t.Contains, htsio.TagContains}, {t.NotContains, htsio.TagNotContains},
		{t.Lt, htsio.TagLt}, {t.Gt, htsio.TagGt},
		{t.Lte, htsio.TagLte}, {t.Gte, htsio.TagGte},
	} {
		for _, spec := range *pair.specs {
			f, err := htsio.ParseTagFilter(spec, pair.op)
			if err != nil {
				return nil, err
			}
			filters = append(filters, f)
		}
	}
	return filters, nil
}

// samReaderFlags holds the common filtering flags shared across SAM commands.
type samReaderFlags struct {
	flagRequired int
	flagExclude  int
	minMapQ      int
	region       string
	ref          string
	tags         tagFilterFlags
}

// register adds the common filtering flags to the command.
func (f *samReaderFlags) register(cmd *cobra.Command) {
	cmd.Flags().IntVar(&f.flagRequired, "flag-required", 0, "Keep reads with all of these flag bits set")
	cmd.Flags().IntVar(&f.flagExclude, "flag-exclude", 0, "Exclude reads with any of these flag bits set")
	cmd.Flags().IntVar(&f.minMapQ, "min-mapq", 0, "Keep reads with mapping quality at or above this value")
	cmd.Flags().StringVar(&f.region, "region", "", "Genomic region (chrom:start-end)")
	cmd.Flags().StringVar(&f.ref, "ref", "", "Filter by reference name")
	f.tags.register(cmd)
}

// buildReaderOpts creates a SamReaderOpts from the parsed flag values.
func (f *samReaderFlags) buildReaderOpts() (*htsio.SamReaderOpts, error) {
	tagFilters, err := f.tags.parse()
	if err != nil {
		return nil, err
	}

	opts := htsio.NewSamReaderOpts()
	if f.flagRequired != 0 {
		opts.FlagRequired(f.flagRequired)
	}
	if f.flagExclude != 0 {
		opts.FlagFilter(f.flagExclude)
	}
	if f.minMapQ > 0 {
		opts.MinMapQ(f.minMapQ)
	}
	if f.region != "" {
		opts.Region(f.region)
	} else if f.ref != "" {
		opts.Region(f.ref)
	}
	for _, tf := range tagFilters {
		opts.AddTagFilter(tf)
	}
	return opts, nil
}
