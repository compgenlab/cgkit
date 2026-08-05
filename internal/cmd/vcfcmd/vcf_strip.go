package vcfcmd

import (
	"fmt"
	"os"
	"strings"

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
	vcfStripForceEnd   bool
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
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
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
		infoSet, err := newStripSet("info", removeInfo, vcfStripKeepInfo)
		if err != nil {
			return err
		}
		formatSet, err := newStripSet("format", removeFormat, vcfStripKeepFormat)
		if err != nil {
			return err
		}
		filterSet, err := newStripSet("filter", removeFilter, vcfStripKeepFilter)
		if err != nil {
			return err
		}
		sampleSet, err := newStripSet("sample", removeSample, vcfStripKeepSample)
		if err != nil {
			return err
		}

		var keptIdx []int
		var keepEnd bool
		return runVcfStream(cmd, vcfStream{
			name: "vcf-strip",
			in:   args[0],
			out:  vcfStripOutput,
			header: func(header *vcf.VcfHeader) error {
				// A gVCF's reference blocks declare their extent with INFO/END. Removing it
				// leaves a file that still parses but silently claims one base of coverage
				// where it claimed thousands, so rescue it unless told otherwise.
				gvcf := isGvcfHeader(header)
				keepEnd = gvcf && infoSet.strips("END") && !vcfStripForceEnd
				if gvcf && infoSet.strips("END") {
					warnStrippingGvcfEnd(cmd, keepEnd)
				}

				// Header: drop stripped INFO/FORMAT/FILTER defs.
				for _, id := range append([]string(nil), header.InfoIDs()...) {
					if keepEnd && id == "END" {
						continue
					}
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
				var keptNames []string
				for i, s := range header.Samples() {
					if !sampleSet.strips(s) {
						keptIdx = append(keptIdx, i)
						keptNames = append(keptNames, s)
					}
				}
				header.SetSamples(keptNames)
				return nil
			},
			record: func(rec *vcf.VcfRecord) (bool, error) {
				if err := stripRecord(rec, infoSet, formatSet, filterSet, keptIdx, dbsnp, keepEnd, vcfStripForceEnd); err != nil {
					return false, err
				}
				switch {
				case vcfStripPassing && rec.IsFiltered():
					return false, nil
				case vcfStripOnlySNVs && rec.IsIndel():
					return false, nil
				case vcfStripOnlyIndels && !rec.IsIndel():
					return false, nil
				}
				return true, nil
			},
		})
	},
}

func stripRecord(rec *vcf.VcfRecord, infoSet, formatSet, filterSet stripSet, keptIdx []int, dbsnp bool, keepEnd, forceEnd bool) error {
	// Last line of defence for a gVCF whose header declared nothing that
	// isGvcfHeader recognises. By now the output header has already been written,
	// so END cannot be reinstated -- and continuing would emit records carrying an
	// undeclared END, or none at all. Fail instead of writing that file.
	if !keepEnd && !forceEnd && infoSet.strips("END") && rec.IsRefBlock() {
		if _, ok := rec.Info().Get("END"); ok {
			return fmt.Errorf("%s:%d is a reference block with INFO/END, but END is being stripped\n"+
				"       and the header did not identify this file as a gVCF.\n"+
				"       Removing END would silently reduce every block to one base.\n"+
				"       Use --keep-info END to retain it, or --force-end to remove it anyway",
				rec.Chrom, rec.Pos)
		}
	}
	// INFO keys.
	for _, k := range append([]string(nil), rec.Info().Keys()...) {
		if keepEnd && k == "END" {
			continue
		}
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

// newStripSet resolves a remove/keep flag pair, expanding any "@path" value
// into that file's lines. flag names the removing flag, for error messages.
func newStripSet(flag string, remove, keep []string) (stripSet, error) {
	r, err := expandStripValues(flag, remove)
	if err != nil {
		return stripSet{}, err
	}
	k, err := expandStripValues("keep-"+flag, keep)
	if err != nil {
		return stripSet{}, err
	}
	return stripSet{remove: r, keep: k}, nil
}

// expandStripValues resolves each --info/--format/--filter/--sample value,
// reading "@path" as a file of names and taking anything else literally.
//
// The "@" is required, though the help has always documented it. The code used
// to probe with os.Stat instead, so the documented syntax did not work at all
// -- "@fields.txt" was taken as a field literally named "@fields.txt" -- while
// an undocumented bare filename did. Worse, that probe made behaviour depend on
// the working directory: these values are short tokens like AC or DP, so a file
// with such a name sitting nearby silently turned a field removal into a list
// read, with nothing said either way.
//
// A bare value naming an existing file is therefore an error rather than a
// silent literal, so anyone relying on the old form is told what to type.
func expandStripValues(flag string, vals []string) ([]string, error) {
	if len(vals) == 0 {
		return nil, nil
	}
	var out []string
	for _, v := range vals {
		if name, ok := strings.CutPrefix(v, "@"); ok {
			lines, err := readLines(name)
			if err != nil {
				return nil, fmt.Errorf("--%s %s: %w", flag, v, err)
			}
			out = append(out, lines...)
			continue
		}
		if fi, err := os.Stat(v); err == nil && !fi.IsDir() {
			return nil, fmt.Errorf("--%s %q names an existing file; write --%s @%s to "+
				"read the names from it, or move the file if you meant the literal name",
				flag, v, flag, v)
		}
		out = append(out, v)
	}
	return out, nil
}

func init() {
	f := vcfStripCmd.Flags()
	addVcfOutputFlags(vcfStripCmd, &vcfStripOutput)
	f.BoolVar(&vcfStripAll, "all", false, "Remove all annotations and samples")
	f.BoolVar(&vcfStripDBSNP, "dbsnp", false, "Remove the ID column")
	f.BoolVar(&vcfStripPassing, "passing", false, "Only output passing variants (post-strip)")
	f.BoolVar(&vcfStripOnlySNVs, "only-snvs", false, "Only output SNVs")
	f.BoolVar(&vcfStripOnlyIndels, "only-indels", false, "Only output indels")
	f.BoolVar(&vcfStripForceEnd, "force-end", false,
		"Allow removing INFO/END from a gVCF (by default it is kept, since removing it collapses every reference block to one base)")
	f.StringArrayVar(&vcfStripInfo, "info", nil, "Remove these INFO fields (glob or @file; repeatable)")
	f.StringArrayVar(&vcfStripFormat, "format", nil, "Remove these FORMAT fields (glob or @file; repeatable)")
	f.StringArrayVar(&vcfStripFilter, "filter", nil, "Remove these FILTER codes (glob or @file; repeatable)")
	f.StringArrayVar(&vcfStripSample, "sample", nil, "Remove these samples (glob or @file; repeatable)")
	f.StringArrayVar(&vcfStripKeepInfo, "keep-info", nil, "Keep these INFO fields (rescue from removal)")
	f.StringArrayVar(&vcfStripKeepFormat, "keep-format", nil, "Keep these FORMAT fields (rescue from removal)")
	f.StringArrayVar(&vcfStripKeepFilter, "keep-filter", nil, "Keep these FILTER codes (rescue from removal)")
	f.StringArrayVar(&vcfStripKeepSample, "keep-sample", nil, "Keep these samples (rescue from removal)")
}
