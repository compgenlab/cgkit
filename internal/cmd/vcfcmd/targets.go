package vcfcmd

import (
	"bufio"
	"context"
	"fmt"
	"io"
	"math"
	"os"
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/iosource"
	"github.com/compgenlab/cghts/varstore"
	"github.com/compgenlab/cghts/vcf"
	"github.com/compgenlab/cgkit/internal/locator"
)

// Target parsing for vcf-varquery's --variant, which is deliberately the ONLY
// site selector: it takes inline loci, inline regions, or a file of either.
//
//	chr1                     a whole contig
//	chr1:1000-2000           a region
//	chr1:1000                any variant at that position
//	chr1:1000:A:T            one exact variant
//	panel.vcf|.bed|.txt      a file, format detected from its content
//
// A value is a file when it carries a URI scheme, or when os.Stat says one
// exists by that name; otherwise it is an inline selector. That order matters
// for error quality: a mistyped locus is not a file, so it still gets a locus
// error rather than "no such file".
//
// The scheme test comes first because a URL cannot be stat'd, and a stat-first
// rule would quietly demote one to an inline selector -- which then parses as a
// contig name and matches nothing, reporting an empty result instead of an
// error. Only values containing "://" past the first character change meaning,
// and nothing in the selector grammar can look like that: a mistyped locus like
// "chr1:100:A" or "chr1:1O0-2000" takes exactly the path it always did, and a
// contig name carrying colons ("HLA-A*01:01:01:01") is untouched.
//
// The site-list grammar is the same one vcf-gtcount's --sites accepts -- whitespace
// separated "chrom pos [ref alt]", '#' comments, blank lines skipped -- so one file
// works with both commands. TestSiteListWorksForBothCommands pins that; the two
// parsers are separate code and the shared FORMAT is what must not drift.

// targetFormat labels where targets came from, for the verbose report. Detection
// is worth reporting because getting it wrong is silent: BED is 0-based
// half-open while site lists and VCF are 1-based, so a misread shifts every
// coordinate by one rather than failing.
type targetFormat string

const (
	targetInline targetFormat = "inline"
	targetVCF    targetFormat = "vcf"
	targetBED    targetFormat = "bed"
	targetList   targetFormat = "site list"
)

// targetRef points at one entry of Loci or Spans. It records the order targets
// were given in, which matters to a consumer that looks each one up in turn and
// reports results in that order -- vcf-gtcount does.
type targetRef struct {
	locus bool
	i     int
}

// targetSet is what --variant resolved to.
type targetSet struct {
	Loci   []varstore.Locus
	Spans  []varstore.Span
	Order  []targetRef          // Loci and Spans interleaved, as given
	Counts map[targetFormat]int // targets contributed by each source kind
	Files  map[string]targetFormat
}

// addLocus and addSpan are the only ways to add a target, so Order cannot drift
// out of step with Loci and Spans.
func (t *targetSet) addLocus(l varstore.Locus) {
	t.Order = append(t.Order, targetRef{locus: true, i: len(t.Loci)})
	t.Loci = append(t.Loci, l)
}

func (t *targetSet) addSpan(sp varstore.Span) {
	t.Order = append(t.Order, targetRef{i: len(t.Spans)})
	t.Spans = append(t.Spans, sp)
}

// parseTargets resolves every --variant value into query selectors.
func parseTargets(ctx context.Context, vals []string) (*targetSet, error) {
	t := &targetSet{Counts: map[targetFormat]int{}, Files: map[string]targetFormat{}}
	for _, v := range vals {
		if isTargetFile(v) {
			if err := t.addFile(ctx, v); err != nil {
				return nil, err
			}
			continue
		}
		if err := t.addInline(v); err != nil {
			return nil, err
		}
	}
	return t, nil
}

// isTargetFile reports whether a --variant or --sample value names a file
// rather than an inline selector. See the package comment above for why the
// scheme test has to precede the stat.
func isTargetFile(v string) bool {
	if locator.IsRemote(v) {
		return true
	}
	st, err := os.Stat(v)
	return err == nil && !st.IsDir()
}

// addInline parses one command-line selector.
func (t *targetSet) addInline(s string) error {
	l, sp, err := parseSelector(s)
	if err != nil {
		return err
	}
	if l != nil {
		t.addLocus(*l)
	}
	if sp != nil {
		t.addSpan(*sp)
	}
	t.Counts[targetInline]++
	return nil
}

// parseSelector parses one inline selector -- the same grammar accepted on the
// command line and on a single-field line of a target file, so a file can just be
// a list of the tokens you would otherwise type.
func parseSelector(s string) (*varstore.Locus, *varstore.Span, error) {
	parts := strings.Split(s, ":")

	// An exact locus: four colon-fields, a numeric position, and REF/ALT that are
	// not numbers. That last test is what lets a contig name carry colons. GRCh38's
	// ALT contigs are named like HLA-A*01:01:01:01, which also splits into four
	// fields -- but its last two are numeric, where a REF or ALT allele never is.
	// Without the test that contig parses as chrom=HLA-A*01, pos=1, ref=01, alt=01:
	// accepted, wrong, and silent.
	if len(parts) == 4 && isNum(parts[1]) && !isNum(parts[2]) && !isNum(parts[3]) {
		l, err := varstore.ParseLocus(s)
		if err != nil {
			return nil, nil, err
		}
		return &l, nil, nil
	}

	// A region, or a single position.
	if len(parts) == 2 {
		if strings.Contains(parts[1], "-") {
			sp, err := varstore.ParseSpan(s)
			if err != nil {
				return nil, nil, fmt.Errorf("invalid target %q: %w", s, err)
			}
			return nil, sp, nil
		}
		if isNum(parts[1]) {
			// chrom:pos means any variant at that position, which is a one-base span
			// -- a Locus cannot express it, since it requires ref and alt.
			pos, err := strconv.Atoi(parts[1])
			if err != nil || pos < 1 {
				return nil, nil, fmt.Errorf("invalid target %q: bad position", s)
			}
			return nil, &varstore.Span{
				Chrom: parts[0], Start: int32(pos - 1), End: int32(pos),
			}, nil
		}
	}

	// Anything else is a contig name, taken whole -- colons included. That covers
	// the shorter ALT spellings too (HLA-B*07:02:01 is three fields, so it is
	// neither a locus nor a region). A typo like "chr1:100:A" also lands here and
	// selects nothing, which is why a selector that matched no rows is reported.
	//
	// A locus ON such a contig cannot be written inline at all, since the grammar
	// counts colons; a target file's columnar form takes it as separate fields.
	return nil, &varstore.Span{Chrom: s, Start: 0, End: math.MaxInt32}, nil
}

// isNum reports whether a field is an integer, which is how a REF/ALT allele is
// told from a numeric component of a contig name or a coordinate.
func isNum(s string) bool {
	_, err := strconv.Atoi(s)
	return err == nil
}

// addFile detects a target file's format and reads it.
func (t *targetSet) addFile(ctx context.Context, path string) error {
	format, err := detectTargetFile(ctx, path)
	if err != nil {
		return err
	}
	t.Files[path] = format
	before := len(t.Loci) + len(t.Spans)
	switch format {
	case targetVCF:
		err = t.addVcfFile(ctx, path)
	case targetBED:
		err = t.addBedFile(ctx, path)
	default:
		err = t.addSiteList(ctx, path)
	}
	if err != nil {
		return err
	}
	n := len(t.Loci) + len(t.Spans) - before
	if n == 0 {
		return fmt.Errorf("%s: read as %s but contained no targets", path, format)
	}
	t.Counts[format] += n
	return nil
}

// detectTargetFile decides what a target file is from its content.
//
// A VCF announces itself. Otherwise the third column decides: a BED's third
// column is an end coordinate and always numeric, where a site list's is a REF
// allele and never is. Fewer than three columns can only be a site list, since a
// BED interval needs three.
func detectTargetFile(ctx context.Context, path string) (targetFormat, error) {
	f, err := iosource.OpenReader(ctx, path)
	if err != nil {
		return "", err
	}
	defer f.Close()

	sc := bufio.NewScanner(f)
	first := true
	for sc.Scan() {
		text := strings.TrimSpace(sc.Text())
		// Checked before comments are skipped, since the marker is itself a '#' line.
		if first && strings.HasPrefix(text, "##fileformat=VCF") {
			return targetVCF, nil
		}
		first = false
		if text == "" || strings.HasPrefix(text, "#") ||
			strings.HasPrefix(text, "track") || strings.HasPrefix(text, "browser") {
			continue
		}
		fields := strings.Fields(text)
		if len(fields) >= 3 {
			if _, err := strconv.Atoi(fields[2]); err == nil {
				return targetBED, nil
			}
		}
		return targetList, nil
	}
	if err := sc.Err(); err != nil {
		return "", err
	}
	return "", fmt.Errorf("%s: no targets found; expected a VCF, a BED, or "+
		"whitespace-separated 'chrom pos [ref alt]' lines", path)
}

// addVcfFile takes one target per ALT allele, matching how a store splits a
// multiallelic record into one row per alternate.
func (t *targetSet) addVcfFile(ctx context.Context, path string) error {
	r, err := vcf.OpenVcfFile(ctx, path)
	if err != nil {
		return err
	}
	defer r.Close()
	if _, err := r.Header(); err != nil {
		return err
	}
	for {
		rec, err := r.NextRecord()
		if err == io.EOF {
			return nil
		}
		if err != nil {
			return err
		}
		loci := varstore.RecordLoci(rec)
		if len(loci) == 0 {
			// ALT is '.', so there is no alternate to ask about. Counted rather than
			// dropped in silence: losing panel rows changes a score's denominator.
			t.Counts["skipped (no ALT)"]++
			continue
		}
		for _, l := range loci {
			t.addLocus(l)
		}
	}
}

// addBedFile reads 0-based half-open intervals, which is what Span already is.
func (t *targetSet) addBedFile(ctx context.Context, path string) error {
	return eachTargetLine(ctx, path, func(fields []string, where string) error {
		if len(fields) < 3 {
			return fmt.Errorf("%s: expected 'chrom start end'", where)
		}
		start, err1 := strconv.Atoi(fields[1])
		end, err2 := strconv.Atoi(fields[2])
		if err1 != nil || err2 != nil {
			return fmt.Errorf("%s: expected numeric start and end", where)
		}
		if end < start {
			return fmt.Errorf("%s: end %d is before start %d", where, end, start)
		}
		t.addSpan(varstore.Span{Chrom: fields[0], Start: int32(start), End: int32(end)})
		return nil
	})
}

// addSiteList reads 1-based "chrom pos [ref alt]" lines. Without ref and alt the
// line names a position rather than a variant, which becomes a one-base span.
func (t *targetSet) addSiteList(ctx context.Context, path string) error {
	return eachTargetLine(ctx, path, func(fields []string, where string) error {
		// One field is a colon-delimited selector -- the same grammar as the command
		// line -- so a file may simply list the tokens you would otherwise type, and
		// may mix them with columnar lines.
		if len(fields) == 1 {
			l, sp, err := parseSelector(fields[0])
			if err != nil {
				return fmt.Errorf("%s: %w", where, err)
			}
			if l != nil {
				t.addLocus(*l)
			}
			if sp != nil {
				t.addSpan(*sp)
			}
			return nil
		}
		if len(fields) < 2 {
			return fmt.Errorf("%s: expected 'chrom pos [ref alt]'", where)
		}
		pos, err := strconv.Atoi(fields[1])
		if err != nil || pos < 1 {
			return fmt.Errorf("%s: expected a 1-based position", where)
		}
		if len(fields) >= 4 {
			t.addLocus(varstore.Locus{
				Chrom: fields[0], Pos: int32(pos), Ref: fields[2], Alt: fields[3],
			})
			return nil
		}
		t.addSpan(varstore.Span{Chrom: fields[0], Start: int32(pos - 1), End: int32(pos)})
		return nil
	})
}

// eachTargetLine walks a whitespace-separated file, skipping blanks, '#' comments
// and BED track lines. Any whitespace separates, not only tabs.
func eachTargetLine(ctx context.Context, path string, fn func(fields []string, where string) error) error {
	f, err := iosource.OpenReader(ctx, path)
	if err != nil {
		return err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	line := 0
	for sc.Scan() {
		line++
		text := strings.TrimSpace(sc.Text())
		if text == "" || strings.HasPrefix(text, "#") ||
			strings.HasPrefix(text, "track") || strings.HasPrefix(text, "browser") {
			continue
		}
		if err := fn(strings.Fields(text), fmt.Sprintf("%s:%d", path, line)); err != nil {
			return err
		}
	}
	return sc.Err()
}

// report describes what --variant resolved to, on stderr.
//
// Worth printing because detection is silent when wrong: a BED read as a site
// list would shift every coordinate by one and still produce plausible output.
func (t *targetSet) report(w io.Writer) {
	for path, format := range t.Files {
		fmt.Fprintf(w, "targets  %s read as %s\n", path, format)
	}
	if n := t.Counts[targetInline]; n > 0 {
		fmt.Fprintf(w, "targets  %d given inline\n", n)
	}
	for _, f := range []targetFormat{targetVCF, targetBED, targetList} {
		if n := t.Counts[f]; n > 0 {
			fmt.Fprintf(w, "targets  %d from %s\n", n, f)
		}
	}
	if n := t.Counts["skipped (no ALT)"]; n > 0 {
		fmt.Fprintf(w, "targets  %d record(s) skipped: ALT is '.', no alternate to ask about\n", n)
	}
	fmt.Fprintf(w, "targets  %d locus/loci, %d span(s)\n", len(t.Loci), len(t.Spans))
}

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
		if !isTargetFile(v) {
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

	err = eachTargetLine(ctx, path, func(fields []string, _ string) error {
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
