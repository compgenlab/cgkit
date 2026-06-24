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
	vcfAnnotateOutput      string
	vcfAnnotatePassing     bool
	vcfAnnotateAltChrom    string
	vcfAnnotateAltPos      string
	vcfAnnotateEndPos      string
	vcfAnnotateAutoID      bool
	vcfAnnotateTags        []string
	vcfAnnotateIndel       bool
	vcfAnnotateTsTv        bool
	vcfAnnotateDosage      bool
	vcfAnnotateVarDist     bool
	vcfAnnotateVAF         bool
	vcfAnnotateMinorStrand bool
	vcfAnnotateFisherSB    bool
	vcfAnnotateCopyLR      string
	vcfAnnotateBed         []string
	vcfAnnotateBedFlag     []string
	vcfAnnotateFormatBed   []string
	vcfAnnotateTab         []string
	vcfAnnotateFormatTab   []string
)

var vcfAnnotateCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-annotate <input.vcf>",
	Short:       "Annotate a VCF file by adding INFO/FORMAT fields",
	Long: `Annotate a VCF file. Each flag adds an annotator that writes one or more
INFO/FORMAT fields onto every record. Annotators are applied in a fixed order
and a matching ##INFO/##FORMAT header line is added for each new field.

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

// buildAnnotatePipeline constructs the annotator pipeline from the flags, in a
// fixed, documented order (cobra cannot recover cross-flag CLI order; the order
// only affects header-line order since the annotations are independent).
func buildAnnotatePipeline() (*annotate.Pipeline, error) {
	p := annotate.NewPipeline()
	add := func(a annotate.Annotator) {
		applyAltCoords(a)
		p.Add(a)
	}

	if vcfAnnotateAutoID {
		add(annotate.NewAutoID())
	}
	for _, tag := range vcfAnnotateTags {
		if i := strings.IndexByte(tag, ':'); i >= 0 {
			add(annotate.NewConstantTag(tag[:i], tag[i+1:]))
		} else {
			add(annotate.NewConstantFlag(tag))
		}
	}
	if vcfAnnotateIndel {
		add(annotate.NewIndel())
	}
	if vcfAnnotateTsTv {
		add(annotate.NewTsTv())
	}
	if vcfAnnotateDosage {
		add(annotate.NewDosage())
	}
	if vcfAnnotateVAF {
		add(annotate.NewVAF())
	}
	if vcfAnnotateMinorStrand {
		add(annotate.NewMinorStrand())
	}
	if vcfAnnotateFisherSB {
		add(annotate.NewFisherSB())
	}
	if vcfAnnotateCopyLR != "" {
		a, err := parseCopyLogRatio(vcfAnnotateCopyLR)
		if err != nil {
			return nil, err
		}
		add(a)
	}

	addTabix := func(o annotate.TabixOptions) error {
		a, err := annotate.NewTabixAnnotator(o)
		if err != nil {
			return err
		}
		applyAltCoords(a)
		p.Add(a)
		return nil
	}
	type tabixSpec struct {
		args  []string
		parse func(string) (annotate.TabixOptions, error)
	}
	for _, spec := range []tabixSpec{
		{vcfAnnotateBed, func(s string) (annotate.TabixOptions, error) { return parseBedArg(s, false) }},
		{vcfAnnotateBedFlag, func(s string) (annotate.TabixOptions, error) { return parseBedArg(s, true) }},
		{vcfAnnotateFormatBed, parseFormatBedArg},
		{vcfAnnotateTab, parseTabArg},
		{vcfAnnotateFormatTab, parseFormatTabArg},
	} {
		for _, arg := range spec.args {
			o, err := spec.parse(arg)
			if err != nil {
				return nil, err
			}
			if err := addTabix(o); err != nil {
				return nil, err
			}
		}
	}

	if vcfAnnotateVarDist {
		s := annotate.NewVariantDistance()
		applyAltCoords(s)
		p.AddStream(s)
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
	f.BoolVar(&vcfAnnotateAutoID, "auto-id", false, "Set the ID to chrom_pos_ref_alt")
	f.StringArrayVar(&vcfAnnotateTags, "tag", nil, "Add a constant INFO annotation: KEY or KEY:VALUE (repeatable)")
	f.BoolVar(&vcfAnnotateIndel, "indel", false, "Add INSERT/DELETE flags and lengths")
	f.BoolVar(&vcfAnnotateTsTv, "tstv", false, "Add TS/TV annotation (CG_TSTV)")
	f.BoolVar(&vcfAnnotateDosage, "dosage", false, "Add per-sample dosage from GT (CG_DS)")
	f.BoolVar(&vcfAnnotateVarDist, "vardist", false, "Add distance to nearest variant (CG_VARDIST)")
	f.BoolVar(&vcfAnnotateVAF, "vaf", false, "Add variant allele frequency (CG_VAF, requires SAC)")
	f.BoolVar(&vcfAnnotateMinorStrand, "minor-strand", false, "Add minor strand percentage (CG_SBPCT, requires SAC)")
	f.BoolVar(&vcfAnnotateFisherSB, "fisher-sb", false, "Add Fisher strand bias (CG_FSB, requires SAC)")
	f.StringVar(&vcfAnnotateCopyLR, "copy-logratio", "", "Add copy-number log2 ratio: SOMATIC:GERMLINE[:somatic-total:germline-total] (CG_CNLR, requires AD)")
	f.StringArrayVar(&vcfAnnotateBed, "bed", nil, "Annotate INFO from a tabix-indexed BED4 name column: NAME:FILE (',n' on NAME for numeric; repeatable)")
	f.StringArrayVar(&vcfAnnotateBedFlag, "bed-flag", nil, "Flag variants within a tabix-indexed BED region: NAME:FILE (repeatable)")
	f.StringArrayVar(&vcfAnnotateFormatBed, "format-bed", nil, "Annotate a sample FORMAT field from a BED4 name column: KEY:SAMPLE:FILE (repeatable)")
	f.StringArrayVar(&vcfAnnotateTab, "tab", nil, "Annotate INFO from a tabix file: NAME:FILE{,col,n,alt=C,ref=C,collapse,first,max,extend=N} (col/alt/ref may be a 1-based number or a header column name when the file has a skipped header line; repeatable)")
	f.StringArrayVar(&vcfAnnotateFormatTab, "format-tab", nil, "Annotate a sample FORMAT field from a tabix file: NAME:SAMPLE:FILE,col{,...} (see --tab; repeatable)")
}
