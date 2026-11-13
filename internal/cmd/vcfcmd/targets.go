package vcfcmd

import (
	"bufio"
	"context"
	"fmt"
	"io"
	"strings"

	"github.com/compgenlab/cghts/iosource"
	"github.com/compgenlab/cghts/varstore"
	"github.com/compgenlab/cghts/vcf"
)

// Site selection (--variant) is parsed by cghts: varstore.ParseTargets in
// varstore/targets.go. What stays here is the sample axis.

// Sample-axis parsing. --sample takes a name or a file of names, mirroring how
// --variant takes a locus or a file of targets.

// sampleSet is what --sample resolved to.
type sampleSet struct {
	Names []string
	Files map[string]string // path -> how it was read, for the verbose report
}

// parseSampleArgs resolves --sample values. A value is a file when one exists by
// that name, and a sample name otherwise -- the same precedence --variant uses.
//
// A VCF is recognised and its header roster taken, which is both useful (subset to
// the samples another callset has) and defensive: read as a name list, a VCF would
// silently yield thousands of bogus "samples" from its data lines.
//
// Names are deduplicated, first occurrence winning, because a repeated sample would
// otherwise appear twice as a column in --format vcf output.
func parseSampleArgs(ctx context.Context, vals []string) (*sampleSet, error) {
	set := &sampleSet{Files: map[string]string{}}
	seen := map[string]bool{}
	add := func(name string) {
		if name != "" && !seen[name] {
			seen[name] = true
			set.Names = append(set.Names, name)
		}
	}

	for _, v := range vals {
		if !varstore.IsTargetFile(v) {
			add(v)
			continue
		}
		names, how, err := readSampleFile(ctx, v)
		if err != nil {
			return nil, err
		}
		if len(names) == 0 {
			return nil, fmt.Errorf("%s: read as %s but named no samples", v, how)
		}
		set.Files[v] = how
		for _, n := range names {
			add(n)
		}
	}
	return set, nil
}

// readSampleFile reads a sample roster: a VCF's header, or whitespace-separated
// names with '#' comments.
func readSampleFile(ctx context.Context, path string) (names []string, how string, err error) {
	f, err := iosource.OpenReader(ctx, path)
	if err != nil {
		return nil, "", err
	}
	sc := bufio.NewScanner(f)
	isVcf := sc.Scan() && strings.HasPrefix(strings.TrimSpace(sc.Text()), "##fileformat=VCF")
	f.Close()

	if isVcf {
		r, err := vcf.OpenVcfFile(ctx, path)
		if err != nil {
			return nil, "", err
		}
		defer r.Close()
		h, err := r.Header()
		if err != nil {
			return nil, "", err
		}
		return h.Samples(), "a VCF header", nil
	}

	err = eachListLine(ctx, path, func(fields []string) error {
		names = append(names, fields...)
		return nil
	})
	return names, "a sample list", err
}

// checkSamples reports requested samples the source does not have.
//
// Worth erroring on rather than ignoring. A store answers an unknown sample with
// silence, so a typo -- in a name or in the path of a sample file, which is then
// taken as a name -- would look exactly like a subject that genuinely carries
// nothing. The VCF backend used to fail with "sample not found" before the query
// API stopped distinguishing them; this restores that at the CLI.
func checkSamples(store varstore.Store, want []string) error {
	if len(want) == 0 {
		return nil
	}
	roster, err := store.Samples()
	if err != nil {
		return err
	}
	have := make(map[string]bool, len(roster))
	for _, n := range roster {
		have[n] = true
	}
	var missing []string
	for _, n := range want {
		if !have[n] {
			missing = append(missing, n)
		}
	}
	if len(missing) == 0 {
		return nil
	}
	return fmt.Errorf("not in this source: %s\n"+
		"       it has %d sample(s)%s", strings.Join(missing, ", "), len(roster),
		firstFew(roster))
}

// firstFew names a few samples, so the error is actionable on a cohort without
// printing thousands of ids.
func firstFew(roster []string) string {
	const max = 6
	if len(roster) == 0 {
		return ""
	}
	shown := roster
	suffix := ""
	if len(shown) > max {
		shown, suffix = shown[:max], ", ..."
	}
	return ": " + strings.Join(shown, ", ") + suffix
}

// report describes what --sample resolved to, on stderr.
func (s *sampleSet) report(w io.Writer) {
	for path, how := range s.Files {
		fmt.Fprintf(w, "samples  %s read as %s\n", path, how)
	}
	if len(s.Names) > 0 {
		fmt.Fprintf(w, "samples  %d requested\n", len(s.Names))
	}
}

// eachListLine walks a whitespace-separated list file, skipping blanks and '#'
// comments -- the same rules varstore.ParseTargets applies to a site list, so a
// sample list and a site list read alike.
func eachListLine(ctx context.Context, path string, fn func(fields []string) error) error {
	f, err := iosource.OpenReader(ctx, path)
	if err != nil {
		return err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		text := strings.TrimSpace(sc.Text())
		if text == "" || strings.HasPrefix(text, "#") {
			continue
		}
		if err := fn(strings.Fields(text)); err != nil {
			return err
		}
	}
	return sc.Err()
}
