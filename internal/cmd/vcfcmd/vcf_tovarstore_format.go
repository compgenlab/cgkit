package vcfcmd

import (
	"fmt"
	"sort"
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/varstore"
	"github.com/compgenlab/cghts/vcf"
)

// --format: capturing per-sample FORMAT fields onto the ALT calls.
//
// The counterpart to --info, and the cardinality fits better: a call row is one
// sample at one ALT, so Number=A maps to it exactly.
//
// AN OPTION, NOT A DEFAULT. calls.parquet is the large table -- roughly a
// hundred rows for every one in the sites catalog -- so a column here is a
// hundred times the cost of the same column there. DP, GQ and AD are stored on
// every call already; this is for the ones that are not.

// resolveFormatFields turns --format selectors into declared fields, reading
// each one's Type and Number from the header.
//
// Same rules as --info: types come from the file rather than the command line,
// a glob may match nothing storable and skips it, and naming a field explicitly
// that cannot be stored is an error rather than a silent omission.
func resolveFormatFields(h *vcf.VcfHeader, selectors []string) ([]varstore.FormatField, []string, error) {
	if len(selectors) == 0 {
		return nil, nil, nil
	}
	seen := map[string]bool{}
	var out []varstore.FormatField
	var skipped []string

	for _, raw := range selectors {
		for _, sel := range strings.Split(raw, ",") {
			sel = strings.TrimSpace(sel)
			if sel == "" {
				continue
			}
			ids := []string{sel}
			glob := strings.ContainsAny(sel, "*?")
			if glob {
				ids = h.MatchFormatIDs(sel)
				if len(ids) == 0 {
					return nil, nil, fmt.Errorf("--format %s: no FORMAT field in the header matches", sel)
				}
			}
			for _, id := range ids {
				if seen[id] {
					continue
				}
				def, ok := h.FormatDef(id)
				if !ok {
					return nil, nil, fmt.Errorf(
						"--format %s: the header declares no ##FORMAT field by that name%s",
						id, nearestFormat(h, id))
				}
				f := varstore.FormatField{
					Name:   id,
					Column: varstore.FormatColumn(id),
					Type:   varstore.InfoType(def.Type),
					Number: def.Number,
				}
				if err := varstore.ValidateFormat([]varstore.FormatField{f}); err != nil {
					// Naming a field is a request; matching one is a side
					// effect. Without the distinction `--format '*'` could never
					// succeed, since every VCF declares GT and most declare PL.
					if glob {
						skipped = append(skipped, fmt.Sprintf("%s (Number=%s)", id, def.Number))
						continue
					}
					return nil, nil, err
				}
				seen[id] = true
				out = append(out, f)
			}
		}
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Name < out[j].Name })
	return out, skipped, nil
}

func nearestFormat(h *vcf.VcfHeader, want string) string {
	for _, id := range h.FormatIDs() {
		if strings.EqualFold(id, want) {
			return fmt.Sprintf(" (it does declare %s -- FORMAT keys are case sensitive)", id)
		}
	}
	return ""
}

// captureFormat reads one sample's values for the captured fields into dst,
// which is reused between calls.
//
// A key the sample does not carry is LEFT OUT rather than zeroed, so the store
// writes a null: "this sample published no VAF here" and "its VAF is 0" are
// different claims.
func captureFormat(
	dst map[string]any, rec *vcf.VcfRecord, sampleIdx int,
	fields []varstore.FormatField, altIdx int,
) {
	clear(dst)
	if len(fields) == 0 {
		return
	}
	attrs, err := rec.Sample(sampleIdx)
	if err != nil {
		return
	}
	for _, f := range fields {
		v, ok := attrs.Get(f.Name)
		if !ok {
			continue
		}
		text := v.String()
		if text == "" || text == "." {
			continue
		}
		if f.Number == "A" {
			parts := strings.Split(text, ",")
			if altIdx >= len(parts) {
				// Fewer values than ALTs: attaching the first to another
				// allele's row would be a number that means something else.
				continue
			}
			text = parts[altIdx]
		}
		if text == "" || text == "." {
			continue
		}
		switch f.Type {
		case varstore.InfoInteger:
			n, err := strconv.ParseInt(text, 10, 32)
			if err != nil {
				continue
			}
			dst[f.Name] = int32(n)
		case varstore.InfoFloat:
			x, err := strconv.ParseFloat(text, 64)
			if err != nil {
				continue
			}
			dst[f.Name] = x
		case varstore.InfoString:
			dst[f.Name] = text
		}
	}
}
