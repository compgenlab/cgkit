package vcfcmd

import (
	"fmt"
	"strings"

	"github.com/spf13/cobra"

	"github.com/compgenlab/cghts/vcf"
)

// isGvcfHeader reports whether a header describes a gVCF -- a file carrying
// reference blocks, not just variants.
//
// Detection has to work from the header alone, because `vcf-strip` writes its
// output header before reading the first record, so a record-level answer arrives
// too late to keep the ##INFO declaration. None of these markers is mandated by the
// spec, so any one is taken as sufficient:
//
//   - a declared <NON_REF> (GATK) or <*> (VCF 4.5) alternate
//   - a ##GVCFBlock... line, which GATK writes to describe its GQ banding
//   - a declared FORMAT/MIN_DP or FORMAT/LEN, both of which only mean anything for
//     a reference block
//
// A gVCF that declares none of these is possible, which is why callers should also
// fail safe when they meet a reference block they were not expecting.
//
// None of these markers is *evidence* of reference blocks, only of a file whose
// header has at some point described them, so this is far too weak a basis for a
// refusal -- see gvcfRefBlockError. A DRAGEN msVCF, which is a joint-genotyped
// cohort callset and the exact input vcf-tovarstore is for, keeps the
// ##ALT=<ID=NON_REF> line it inherited from the gVCFs it was built from. VCF 4.5
// makes the <*> case worse than incidental: the gVCF unspecified allele and the
// ordinary spanning-deletion allele are the same ALT ID, so no header-level test
// can separate them even in principle.
func isGvcfHeader(h *vcf.VcfHeader) bool {
	if h == nil {
		return false
	}
	for _, id := range h.AltIDs() {
		if id == "NON_REF" || id == "*" {
			return true
		}
	}
	for _, id := range h.FormatIDs() {
		if id == "MIN_DP" || id == "LEN" {
			return true
		}
	}
	for _, line := range h.OtherLines() {
		if strings.HasPrefix(line, "##GVCFBlock") {
			return true
		}
	}
	return false
}

// gvcfRefBlockError refuses a conversion that has met an actual reference block.
//
// vcf-tovarstore asks this of records rather than of the header, unlike vcf-strip,
// and the difference is not a preference: the header says only what a file's
// ancestry declared, while a record is the thing the store would get wrong. A
// header-based refusal rejected ordinary cohort VCFs -- any callset carrying a
// vestigial ##ALT=<ID=NON_REF>, or a plain ##ALT=<ID=*> for spanning deletions --
// with the store's three-way wrongness argued about a file that has no blocks in it.
// Reading first costs nothing recoverable, because a failed conversion Discards.
//
// The test is IsRefBlock, so it is the *block* that is refused, not the allele: a
// mixed "G,<NON_REF>" record is a variant record and converts, with the block
// allele masked out the way any non-focal alternate is.
func gvcfRefBlockError(path string, rec *vcf.VcfRecord) error {
	return fmt.Errorf("%s contains gVCF reference blocks (%s:%d %s>%s), and converting one is not supported yet\n"+
		"       They would be stored as if they were variants:\n"+
		"         - <NON_REF> would enter the sites catalog, which is meant to be\n"+
		"           the exact boundary of what the store can answer\n"+
		"         - AC/AN would count a reference-block allele as an allele\n"+
		"         - each block's span would be discarded, keeping only its first base\n"+
		"       Query the gVCF directly with vcf-varquery, which reads blocks as the\n"+
		"       coverage they are.", path, rec.Chrom, rec.Pos, rec.Ref, rec.AltOrig())
}

// warnStrippingGvcfEnd explains what removing INFO/END from a gVCF would do.
//
// Not gated behind -v, and printed in both directions. A reference block's extent
// lives entirely in END: without it the record still parses and still says
// "reference here", but for one base instead of thousands. That is a silent change
// in what the file asserts, which is exactly the case a user needs told about.
func warnStrippingGvcfEnd(cmd *cobra.Command, kept bool) {
	if kept {
		fmt.Fprintf(cmd.ErrOrStderr(),
			"warning: this looks like a gVCF, so INFO/END was kept despite the strip options.\n"+
				"         Removing it would leave every reference block claiming one base of\n"+
				"         coverage instead of its real extent. Pass --force-end to remove it.\n")
		return
	}
	fmt.Fprintf(cmd.ErrOrStderr(),
		"warning: removing INFO/END from a gVCF as requested by --force-end.\n"+
			"         Every reference block in the output will claim one base of coverage\n"+
			"         instead of its real extent.\n")
}
