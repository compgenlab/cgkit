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
}
