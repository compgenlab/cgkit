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
// Detection has to work from the header alone, because the commands that need it
// write their output header before reading the first record, so a record-level
// answer arrives too late to act on. None of these markers is mandated by the
// spec, so any one is taken as sufficient:
//
//   - a declared <NON_REF> (GATK) or <*> (VCF 4.5) alternate
//   - a ##GVCFBlock... line, which GATK writes to describe its GQ banding
//   - a declared FORMAT/MIN_DP or FORMAT/LEN, both of which only mean anything for
//     a reference block
//
// A gVCF that declares none of these is possible, which is why callers should also
// fail safe when they meet a reference block they were not expecting.
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
