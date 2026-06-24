package vcfcmd

import (
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/compgenlab/hts/vcf"
	"github.com/compgenlab/hts/vcf/annotate"
	"github.com/spf13/cobra"
)

var (
	vcfAnnotateOutput   string
	vcfAnnotatePassing  bool
	vcfAnnotateAltChrom string
	vcfAnnotateAltPos   string
	vcfAnnotateEndPos   string

	// vcfAnnotateChain records the annotator flags in command-line order so the
	// pipeline is built in that order. See chainValue.
	vcfAnnotateChain []chainArg
)

// chainArg is one annotator flag occurrence: the flag name (kind) and its
// argument (value is "true" for boolean flags).
type chainArg struct {
	kind  string
	value string
}

// chainValue is a pflag.Value that appends to vcfAnnotateChain when set. Because
// pflag calls Set in command-line order, the chain captures the exact order (and
// interleaving) of the annotator flags.
type chainValue struct {
	kind   string
	isBool bool
}

func (c *chainValue) String() string { return "" }
func (c *chainValue) Set(v string) error {
	vcfAnnotateChain = append(vcfAnnotateChain, chainArg{kind: c.kind, value: v})
	return nil
}
func (c *chainValue) Type() string {
	if c.isBool {
		return "bool"
	}
	return "string"
}

// IsBoolFlag makes pflag treat the flag as a no-argument boolean when isBool.
func (c *chainValue) IsBoolFlag() bool { return c.isBool }

var vcfAnnotateCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-annotate <input.vcf>",
	Short:       "Annotate a VCF file by adding INFO/FORMAT fields",
	Long: `Annotate a VCF file. Each flag adds an annotator that writes one or more
INFO/FORMAT fields onto every record. Annotators are applied in the order the
flags appear on the command line (so a later annotator can use a field an
earlier one added), and a matching ##INFO/##FORMAT header line is added for each
new field.

Self-contained annotators (read only the variant):
  --auto-id        set ID to chrom_pos_ref_alt
  --tag KEY[:VAL]  add a constant INFO flag/value (repeatable)
  --indel          flag insertions/deletions and their lengths
  --tstv           CG_TSTV transition/transversion class
  --dosage         CG_DS per-sample dosage from GT
  --vardist        CG_VARDIST distance to nearest variant (sorted input)

Sample-count annotators (require GATK-style FORMAT fields):
  --vaf            CG_VAF allele frequency (requires SAC)
  --minor-strand   CG_SBPCT minor-strand percentage (requires SAC)
  --fisher-sb      CG_FSB Fisher strand bias, Phred-scaled (requires SAC)
  --copy-logratio SOMATIC:GERMLINE[:somatic-total:germline-total]
                   CG_CNLR copy-number log2 ratio (requires AD)`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}

		pipeline, err := buildAnnotatePipeline()
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
		if err := pipeline.SetupHeaders(header); err != nil {
			return err
		}
		stampVcfProvenance(header, "vcf-annotate")

		var writer *vcf.VcfWriter
		var closeFile func() error
		if vcfAnnotateOutput == "" || vcfAnnotateOutput == "-" {
			writer = vcf.NewVcfWriter(cmd.OutOrStdout())
		} else {
			w, err := vcf.OpenVcfWriter(vcfAnnotateOutput)
			if err != nil {
				return err
			}
			writer = w
			closeFile = w.Close
		}
		if err := writer.WriteHeader(header); err != nil {
			return err
		}

		// Source: stream records, optionally dropping filtered ones.
		source := func() (*vcf.VcfRecord, error) {
			for {
				rec, err := reader.NextRecord()
				if err != nil {
					return nil, err
				}
				if vcfAnnotatePassing && rec.IsFiltered() {
					continue
				}
				return rec, nil
			}
		}
		next := pipeline.Build(source)

		for {
			rec, err := next()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if err := writer.WriteRecord(rec); err != nil {
				return err
			}
		}
		if err := pipeline.Close(); err != nil {
			return err
		}
		if closeFile != nil {
			return closeFile()
		}
		return writer.Close()
	},
}

// buildAnnotatePipeline constructs the annotator pipeline from the chain of
// flags in command-line order (see chainValue).
func buildAnnotatePipeline() (*annotate.Pipeline, error) {
	p := annotate.NewPipeline()
	add := func(a annotate.Annotator) {
		applyAltCoords(a)
		p.Add(a)
	}
	// addOpened adds an annotator from a constructor that may fail to open a file.
	addOpened := func(a annotate.Annotator, err error) error {
		if err != nil {
			return err
		}
		add(a)
		return nil
	}
	// addTabix builds a tabix annotator from a parsed-options result.
	addTabix := func(o annotate.TabixOptions, perr error) error {
		if perr != nil {
			return perr
		}
		return addOpened(annotate.NewTabixAnnotator(o))
	}

	for _, c := range vcfAnnotateChain {
		var err error
		switch c.kind {
		case "auto-id":
			add(annotate.NewAutoID())
		case "indel":
			add(annotate.NewIndel())
		case "tstv":
			add(annotate.NewTsTv())
		case "dosage":
			add(annotate.NewDosage())
		case "vaf":
			add(annotate.NewVAF())
		case "minor-strand":
			add(annotate.NewMinorStrand())
		case "fisher-sb":
			add(annotate.NewFisherSB())
		case "tag":
			if i := strings.IndexByte(c.value, ':'); i >= 0 {
				add(annotate.NewConstantTag(c.value[:i], c.value[i+1:]))
			} else {
				add(annotate.NewConstantFlag(c.value))
			}
		case "copy-logratio":
			a, e := parseCopyLogRatio(c.value)
			if e != nil {
				return nil, e
			}
			add(a)
		case "bed":
			err = addTabix(parseBedArg(c.value, false))
		case "bed-flag":
			err = addTabix(parseBedArg(c.value, true))
		case "format-bed":
			err = addTabix(parseFormatBedArg(c.value))
		case "tab":
			err = addTabix(parseTabArg(c.value))
		case "format-tab":
			err = addTabix(parseFormatTabArg(c.value))
		case "vcf":
			o, e := parseVcfArg(c.value)
			if e != nil {
				return nil, e
			}
			err = addOpened(annotate.NewVcfAnnotation(o))
		case "vcf-flag":
			o, e := parseVcfFlagArg(c.value)
			if e != nil {
				return nil, e
			}
			err = addOpened(annotate.NewVcfAnnotation(o))
		case "vcf-id":
			err = addOpened(annotate.NewVcfAnnotation(annotate.VcfOptions{Name: "@ID", Filename: c.value}))
		case "in-file":
			o, e := parseInFileArg(c.value)
			if e != nil {
				return nil, e
			}
			err = addOpened(annotate.NewInfoInFile(o))
		case "vardist":
			s := annotate.NewVariantDistance()
			applyAltCoords(s)
			p.AddStream(s)
		}
		if err != nil {
			return nil, err
		}
	}
	return p, nil
}

// applyAltCoords passes the global --alt-chrom/--alt-pos/--end-pos overrides to
// annotators that resolve query coordinates (none of the current group A/B
// annotators do, so this is currently inert; it is wired for the upcoming
// external-file annotators).
func applyAltCoords(a any) {
	c, ok := a.(annotate.CoordAware)
	if !ok {
		return
	}
	if vcfAnnotateAltChrom != "" {
		c.SetAltChrom(vcfAnnotateAltChrom)
	}
	if vcfAnnotateAltPos != "" {
		c.SetAltPos(vcfAnnotateAltPos)
	}
	if vcfAnnotateEndPos != "" {
		c.SetEndPos(vcfAnnotateEndPos)
	}
}

// splitNumericName strips a trailing ",n" (numeric marker) from a NAME token.
func splitNumericName(name string) (string, bool) {
	if strings.HasSuffix(name, ",n") {
		return name[:len(name)-2], true
	}
	return name, false
}

// parseBedArg parses "NAME:FILE" for --bed / --bed-flag. A flag annotation uses
// Col=0; otherwise the BED name column (4). A ",n" suffix on NAME makes the
// value numeric.
func parseBedArg(arg string, flag bool) (annotate.TabixOptions, error) {
	parts := strings.SplitN(arg, ":", 2)
	if len(parts) != 2 {
		return annotate.TabixOptions{}, fmt.Errorf("expected NAME:FILE, got %q", arg)
	}
	name, isNum := splitNumericName(parts[0])
	o := annotate.TabixOptions{Name: name, Filename: parts[1], IsNumber: isNum, Col: 4}
	if flag {
		o.Col = 0
	}
	return o, nil
}

// parseFormatBedArg parses "KEY:SAMPLE:FILE" for --format-bed.
func parseFormatBedArg(arg string) (annotate.TabixOptions, error) {
	parts := strings.SplitN(arg, ":", 3)
	if len(parts) != 3 {
		return annotate.TabixOptions{}, fmt.Errorf("expected KEY:SAMPLE:FILE, got %q", arg)
	}
	name, isNum := splitNumericName(parts[0])
	return annotate.TabixOptions{Name: name, Sample: parts[1], Filename: parts[2], IsNumber: isNum, Col: 4}, nil
}

// parseTabArg parses "NAME:FILE,opt,..." for --tab.
func parseTabArg(arg string) (annotate.TabixOptions, error) {
	parts := strings.SplitN(arg, ":", 2)
	if len(parts) != 2 {
		return annotate.TabixOptions{}, fmt.Errorf("expected NAME:FILE,..., got %q", arg)
	}
	return parseTabOptions(parts[0], "", parts[1])
}

// parseFormatTabArg parses "NAME:SAMPLE:FILE,opt,..." for --format-tab.
func parseFormatTabArg(arg string) (annotate.TabixOptions, error) {
	parts := strings.SplitN(arg, ":", 3)
	if len(parts) != 3 {
		return annotate.TabixOptions{}, fmt.Errorf("expected NAME:SAMPLE:FILE,..., got %q", arg)
	}
	return parseTabOptions(parts[0], parts[1], parts[2])
}

// parseTabOptions parses the comma-separated FILE,opt,... portion of a --tab /
// --format-tab argument (ports VCFAnnotateCmd.setTabix). Column/alt/ref are
// 1-based numbers; names are not yet supported.
func parseTabOptions(name, sample, fileAndOpts string) (annotate.TabixOptions, error) {
	o := annotate.TabixOptions{Name: name, Sample: sample}
	toks := strings.Split(fileAndOpts, ",")
	colSet := false
	for i, t := range toks {
		switch {
		case i == 0:
			o.Filename = t
		case t == "n":
			o.IsNumber = true
		case t == "max":
			o.Max = true
		case t == "collapse":
			o.Collapse = true
		case t == "first":
			o.First = true
		case t == "noheader":
			o.NoHeader = true
		case strings.HasPrefix(t, "extend="):
			n, err := strconv.Atoi(t[len("extend="):])
			if err != nil {
				return o, fmt.Errorf("invalid extend value: %q", t)
			}
			o.Extend = n
		case strings.HasPrefix(t, "alt="):
			v := t[len("alt="):]
			if n, err := strconv.Atoi(v); err == nil {
				o.AltCol = n
			} else {
				o.AltName = v
			}
		case strings.HasPrefix(t, "ref="):
			v := t[len("ref="):]
			if n, err := strconv.Atoi(v); err == nil {
				o.RefCol = n
			} else {
				o.RefName = v
			}
		case !colSet:
			if n, err := strconv.Atoi(t); err == nil {
				o.Col = n
			} else {
				o.ColName = t
			}
			colSet = true
		}
	}
	if o.Filename == "" {
		return o, fmt.Errorf("missing filename in tab annotation")
	}
	n := 0
	for _, b := range []bool{o.Max, o.First, o.Collapse} {
		if b {
			n++
		}
	}
	if n > 1 {
		return o, fmt.Errorf("first, max, and collapse cannot be combined")
	}
	if o.Max && !o.IsNumber {
		return o, fmt.Errorf("max also requires ,n")
	}
	if o.RefCol > 0 && o.AltCol == 0 {
		return o, fmt.Errorf("ref= requires alt=")
	}
	return o, nil
}

// parseVcfArg parses "NAME:FIELD:FILE[:!@$n]" for --vcf.
func parseVcfArg(arg string) (annotate.VcfOptions, error) {
	spl := strings.Split(arg, ":")
	if len(spl) < 3 {
		return annotate.VcfOptions{}, fmt.Errorf("expected NAME:FIELD:FILE[:mods], got %q", arg)
	}
	o := annotate.VcfOptions{Name: spl[0], Field: spl[1], Filename: spl[2]}
	if len(spl) >= 4 {
		applyVcfMods(&o, spl[3])
	}
	return o, nil
}

// parseVcfFlagArg parses "NAME:FILE[:!@$n]" for --vcf-flag.
func parseVcfFlagArg(arg string) (annotate.VcfOptions, error) {
	spl := strings.Split(arg, ":")
	if len(spl) < 2 {
		return annotate.VcfOptions{}, fmt.Errorf("expected NAME:FILE[:mods], got %q", arg)
	}
	o := annotate.VcfOptions{Name: spl[0], Filename: spl[1]}
	if len(spl) >= 3 {
		applyVcfMods(&o, spl[2])
	}
	return o, nil
}

// parseInFileArg parses "FLAGNAME:INFOKEY:FILE{:csv:tabcol=n}" for --in-file.
func parseInFileArg(arg string) (annotate.InfoFileOptions, error) {
	spl := strings.Split(arg, ":")
	if len(spl) < 3 {
		return annotate.InfoFileOptions{}, fmt.Errorf("expected FLAGNAME:INFOKEY:FILE[:opts], got %q", arg)
	}
	o := annotate.InfoFileOptions{FlagName: spl[0], Tag: spl[1], Filename: spl[2]}
	for _, t := range spl[3:] {
		switch {
		case t == "csv" || t == ",":
			o.Delimiter = ","
		case strings.HasPrefix(t, "tabcol="):
			n, err := strconv.Atoi(t[len("tabcol="):])
			if err != nil {
				return o, fmt.Errorf("invalid tabcol value: %q", t)
			}
			o.Col = n
		}
	}
	return o, nil
}

// applyVcfMods sets the !/@/$/n modifier flags from a modifier string.
func applyVcfMods(o *annotate.VcfOptions, mods string) {
	o.Exact = strings.Contains(mods, "!")
	o.Passing = strings.Contains(mods, "@")
	o.Unique = strings.Contains(mods, "$")
	o.NoHeader = strings.Contains(mods, "n")
}

func parseCopyLogRatio(arg string) (*annotate.CopyNumberLogRatio, error) {
	spl := strings.Split(arg, ":")
	switch len(spl) {
	case 2:
		return annotate.NewCopyLogRatio(spl[0], spl[1], -1, -1), nil
	case 4:
		somTotal, err := strconv.ParseInt(spl[2], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("invalid somatic-total in --copy-logratio: %w", err)
		}
		germTotal, err := strconv.ParseInt(spl[3], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("invalid germline-total in --copy-logratio: %w", err)
		}
		return annotate.NewCopyLogRatio(spl[0], spl[1], somTotal, germTotal), nil
	default:
		return nil, fmt.Errorf("unable to parse --copy-logratio: %q (want SOMATIC:GERMLINE[:somatic-total:germline-total])", arg)
	}
}

func init() {
	f := vcfAnnotateCmd.Flags()
	f.StringVarP(&vcfAnnotateOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
	f.BoolVar(&vcfAnnotatePassing, "passing", false, "Only output passing variants")
	f.StringVar(&vcfAnnotateAltChrom, "alt-chrom", "", "Use an INFO field as the chromosome for coordinate-based annotators")
	f.StringVar(&vcfAnnotateAltPos, "alt-pos", "", "Use an INFO field as the position for coordinate-based annotators")
	f.StringVar(&vcfAnnotateEndPos, "end-pos", "", "Use an INFO field as the end position for coordinate-based annotators")

	// Annotator flags are recorded in command-line order via chainValue. A
	// boolean chain flag needs NoOptDefVal so pflag does not try to consume the
	// next argument as its value.
	chainBool := func(name, usage string) {
		f.Var(&chainValue{kind: name, isBool: true}, name, usage)
		f.Lookup(name).NoOptDefVal = "true"
	}
	chainVal := func(name, usage string) { f.Var(&chainValue{kind: name}, name, usage) }
	chainBool("auto-id", "Set the ID to chrom_pos_ref_alt")
	chainVal("tag", "Add a constant INFO annotation: KEY or KEY:VALUE (repeatable)")
	chainBool("indel", "Add INSERT/DELETE flags and lengths")
	chainBool("tstv", "Add TS/TV annotation (CG_TSTV)")
	chainBool("dosage", "Add per-sample dosage from GT (CG_DS)")
	chainBool("vardist", "Add distance to nearest variant (CG_VARDIST)")
	chainBool("vaf", "Add variant allele frequency (CG_VAF, requires SAC)")
	chainBool("minor-strand", "Add minor strand percentage (CG_SBPCT, requires SAC)")
	chainBool("fisher-sb", "Add Fisher strand bias (CG_FSB, requires SAC)")
	chainVal("copy-logratio", "Add copy-number log2 ratio: SOMATIC:GERMLINE[:somatic-total:germline-total] (CG_CNLR, requires AD)")
	chainVal("bed", "Annotate INFO from a tabix-indexed BED4 name column: NAME:FILE (',n' on NAME for numeric; repeatable)")
	chainVal("bed-flag", "Flag variants within a tabix-indexed BED region: NAME:FILE (repeatable)")
	chainVal("format-bed", "Annotate a sample FORMAT field from a BED4 name column: KEY:SAMPLE:FILE (repeatable)")
	chainVal("tab", "Annotate INFO from a tabix file: NAME:FILE{,col,n,alt=C,ref=C,collapse,first,max,extend=N} (col/alt/ref may be a 1-based number or a header column name when the file has a skipped header line; repeatable)")
	chainVal("format-tab", "Annotate a sample FORMAT field from a tabix file: NAME:SAMPLE:FILE,col{,...} (see --tab; repeatable)")
	chainVal("vcf", "Annotate INFO from a tabix-indexed VCF: NAME:FIELD:FILE{:!@$n} (!=exact ref/alt, @=passing only, $=unique, n=no header def; repeatable)")
	chainVal("vcf-flag", "Flag variants present in a tabix-indexed VCF: NAME:FILE{:!@$n} (repeatable)")
	chainVal("vcf-id", "Copy the ID column from a tabix-indexed VCF (exact ref/alt match)")
	chainVal("in-file", "Flag when an INFO value is present in a text file: FLAGNAME:INFOKEY:FILE{:csv:tabcol=n} (csv splits the INFO value; tabcol=n adds that 1-based column's value; repeatable)")
}
