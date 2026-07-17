package vcfcmd

import (
	"fmt"
	"io"
	"os"

	"github.com/compgenlab/cghts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfStripOutput     string
	vcfStripAll        bool
	vcfStripDBSNP      bool
	vcfStripPassing    bool
	vcfStripOnlySNVs   bool
	vcfStripOnlyIndels bool
	vcfStripInfo       []string
	vcfStripFormat     []string
	vcfStripFilter     []string
	vcfStripSample     []string
	vcfStripKeepInfo   []string
	vcfStripKeepFormat []string
	vcfStripKeepFilter []string
	vcfStripKeepSample []string
)

// stripSet holds the resolved remove/keep glob lists for one field kind.
type stripSet struct{ remove, keep []string }

// strips reports whether id should be removed: it matches a remove glob and is
// not rescued by a keep glob ( porting ngsutilsj's VCFHeader strip logic).
func (s stripSet) strips(id string) bool {
	matched := false
	for _, r := range s.remove {
		if globMatch(id, r) {
			matched = true
			for _, k := range s.keep {
				if globMatch(id, k) {
					matched = false
				}
			}
		}
	}
	return matched
}

var vcfStripCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-strip <input.vcf>",
	Short:       "Remove annotation and sample information, keeping VCF format",
	Long: `Remove annotations (FILTER, INFO, FORMAT, samples, dbSNP ID) from a VCF, while
keeping the output in VCF format.

  --all                       remove all annotations and samples
  --info/--format/--filter/--sample VAL    remove these (glob or @file; repeatable)
  --keep-info/--keep-format/--keep-filter/--keep-sample VAL  rescue these from removal
  --dbsnp                     remove the ID column
  --only-snvs / --only-indels output only SNVs / only indels
  --passing                   output only passing variants (post-strip)`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if vcfStripOnlySNVs && vcfStripOnlyIndels {
			return fmt.Errorf("you can't set both --only-snvs and --only-indels")
		}

		removeInfo, removeFormat := vcfStripInfo, vcfStripFormat
		removeFilter, removeSample := vcfStripFilter, vcfStripSample
		dbsnp := vcfStripDBSNP
		if vcfStripAll {
			dbsnp = true
			removeInfo = append(removeInfo, "*")
			removeFormat = append(removeFormat, "*")
			removeFilter = append(removeFilter, "*")
			removeSample = append(removeSample, "*")
		}
		infoSet, err := newStripSet(removeInfo, vcfStripKeepInfo)
		if err != nil {
			return err
		}
		formatSet, err := newStripSet(removeFormat, vcfStripKeepFormat)
		if err != nil {
			return err
		}
		filterSet, err := newStripSet(removeFilter, vcfStripKeepFilter)
		if err != nil {
			return err
		}
		sampleSet, err := newStripSet(removeSample, vcfStripKeepSample)
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

		// Header: drop stripped INFO/FORMAT/FILTER defs.
		for _, id := range append([]string(nil), header.InfoIDs()...) {
			if infoSet.strips(id) {
				header.RemoveInfo(id)
			}
		}
		for _, id := range append([]string(nil), header.FormatIDs()...) {
			if formatSet.strips(id) {
				header.RemoveFormat(id)
			}
		}
		for _, id := range append([]string(nil), header.FilterIDs()...) {
			if filterSet.strips(id) {
				header.RemoveFilter(id)
			}
		}
		// Header: project samples.
		var keptIdx []int
		var keptNames []string
		for i, s := range header.Samples() {
			if !sampleSet.strips(s) {
				keptIdx = append(keptIdx, i)
				keptNames = append(keptNames, s)
			}
		}
		header.SetSamples(keptNames)
		stampVcfProvenance(header, "vcf-strip")

		writer, closeFn, err := openVcfWriter(cmd, vcfStripOutput)
		if err != nil {
			return err
		}
		if err := writer.WriteHeader(header); err != nil {
			return err
		}

		for {
			rec, err := reader.NextRecord()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if err := stripRecord(rec, infoSet, formatSet, filterSet, keptIdx, dbsnp); err != nil {
				return err
			}
			if vcfStripPassing && rec.IsFiltered() {
				continue
			}
			if vcfStripOnlySNVs && rec.IsIndel() {
				continue
			}
			if vcfStripOnlyIndels && !rec.IsIndel() {
				continue
			}
			if err := writer.WriteRecord(rec); err != nil {
				return err
			}
		}
		if closeFn != nil {
			return closeFn()
		}
		return writer.Close()
	},
}

func stripRecord(rec *vcf.VcfRecord, infoSet, formatSet, filterSet stripSet, keptIdx []int, dbsnp bool) error {
	// INFO keys.
	for _, k := range append([]string(nil), rec.Info().Keys()...) {
		if infoSet.strips(k) {
			rec.Info().Remove(k)
		}
	}
	// FILTER codes (RetainFilters preserves the PASS-vs-"." distinction).
	rec.RetainFilters(func(f string) bool { return !filterSet.strips(f) })
	// FORMAT keys: remove from every kept sample, then project columns.
	var dropKeys []string
	for _, k := range rec.FormatKeys() {
		if formatSet.strips(k) {
			dropKeys = append(dropKeys, k)
		}
	}
	for _, idx := range keptIdx {
		s, err := rec.Sample(idx)
		if err != nil {
			return err
		}
		for _, k := range dropKeys {
			s.Remove(k)
		}
	}
	rec.SubsetSamples(keptIdx)
	if dbsnp {
		rec.ClearID()
	}
	rec.MarkDirty()
	return nil
}

// newStripSet resolves a remove/keep flag pair, expanding any value that names
// an existing file into that file's lines.
func newStripSet(remove, keep []string) (stripSet, error) {
	r, err := expandStripValues(remove)
	if err != nil {
		return stripSet{}, err
	}
	k, err := expandStripValues(keep)
	if err != nil {
		return stripSet{}, err
	}
	return stripSet{remove: r, keep: k}, nil
}

func expandStripValues(vals []string) ([]string, error) {
	if len(vals) == 0 {
		return nil, nil
	}
	var out []string
	for _, v := range vals {
		if fi, err := os.Stat(v); err == nil && !fi.IsDir() {
			lines, err := readLines(v)
			if err != nil {
				return nil, err
			}
			out = append(out, lines...)
		} else {
			out = append(out, v)
		}
	}
	return out, nil
}

func init() {
	f := vcfStripCmd.Flags()
	f.StringVarP(&vcfStripOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
	f.BoolVar(&vcfStripAll, "all", false, "Remove all annotations and samples")
	f.BoolVar(&vcfStripDBSNP, "dbsnp", false, "Remove the ID column")
	f.BoolVar(&vcfStripPassing, "passing", false, "Only output passing variants (post-strip)")
	f.BoolVar(&vcfStripOnlySNVs, "only-snvs", false, "Only output SNVs")
	f.BoolVar(&vcfStripOnlyIndels, "only-indels", false, "Only output indels")
	f.StringArrayVar(&vcfStripInfo, "info", nil, "Remove these INFO fields (glob or @file; repeatable)")
	f.StringArrayVar(&vcfStripFormat, "format", nil, "Remove these FORMAT fields (glob or @file; repeatable)")
	f.StringArrayVar(&vcfStripFilter, "filter", nil, "Remove these FILTER codes (glob or @file; repeatable)")
	f.StringArrayVar(&vcfStripSample, "sample", nil, "Remove these samples (glob or @file; repeatable)")
	f.StringArrayVar(&vcfStripKeepInfo, "keep-info", nil, "Keep these INFO fields (rescue from removal)")
	f.StringArrayVar(&vcfStripKeepFormat, "keep-format", nil, "Keep these FORMAT fields (rescue from removal)")
	f.StringArrayVar(&vcfStripKeepFilter, "keep-filter", nil, "Keep these FILTER codes (rescue from removal)")
	f.StringArrayVar(&vcfStripKeepSample, "keep-sample", nil, "Keep these samples (rescue from removal)")
}
