package vcfcmd

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/compgenlab/hts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfGtCountOutput  string
	vcfGtCountSites   string
	vcfGtCountPassing bool
)

var vcfGtCountCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-gtcount <input.vcf.gz> [locus ...]",
	Short:       "Summarize the genotype (GT) distribution across samples at given sites",
	Long: `For each requested variant site, count how the per-sample GT calls are
distributed across the samples in a multi-sample VCF, writing a tab-delimited
table: chrom, pos, ref, alt, then one column per observed genotype class.

The input VCF must be bgzip-compressed and tabix-indexed (.tbi/.csi); each site
is looked up by an indexed query rather than scanning the whole file.

Sites are given as loci on the command line and/or via --sites:
  locus    chrom:pos              match every record at that position
           chrom:pos:ref:alt      also require REF and ALT to match
  --sites  a file of whitespace-separated columns: chrom pos [ref alt]
           ('#' comments and blank lines are ignored)

Genotypes are collapsed so unordered/phased calls land in one class: 0/1, 1/0
and 0|1 all count as 0/1, and 2/1 as 1/2. Missing calls are reported as ./. (and
absent GT fields count as missing).

  --sites FILE   read additional sites from FILE
  --passing      only count records that pass FILTER`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}

		sites, err := collectGtSites(args[1:], vcfGtCountSites)
		if err != nil {
			return err
		}
		if len(sites) == 0 {
			return fmt.Errorf("no sites given; provide loci as arguments or with --sites")
		}

		ir, err := vcf.NewIndexedVcfReader(args[0])
		if err != nil {
			return err
		}
		defer ir.Close()

		var rows []*gtRow
		classes := map[string]bool{}
		for _, site := range sites {
			matched, err := gtCountSite(ir, site, vcfGtCountPassing)
			if err != nil {
				return err
			}
			if len(matched) == 0 {
				// Emit a zero-count row so every requested site appears.
				matched = []*gtRow{{
					chrom: site.chrom, pos: site.pos,
					ref: orDot(site.ref), alt: orDot(site.alt),
					counts: map[string]int{},
				}}
			}
			for _, r := range matched {
				for c := range r.counts {
					classes[c] = true
				}
				rows = append(rows, r)
			}
		}

		cols := make([]string, 0, len(classes))
		for c := range classes {
			cols = append(cols, c)
		}
		sortGtColumns(cols)

		out, closeFn, err := openOutput(cmd, vcfGtCountOutput)
		if err != nil {
			return err
		}
		fmt.Fprintln(out, "## program: "+buildinfo.String())
		fmt.Fprintln(out, "## cmd: "+buildinfo.CommandLine())
		fmt.Fprintln(out, "## vcf-input: "+args[0])
		fmt.Fprintln(out, strings.Join(append([]string{"chrom", "pos", "ref", "alt"}, cols...), "\t"))

		for _, r := range rows {
			row := []string{r.chrom, strconv.Itoa(r.pos), r.ref, r.alt}
			for _, c := range cols {
				row = append(row, strconv.Itoa(r.counts[c]))
			}
			fmt.Fprintln(out, strings.Join(row, "\t"))
		}

		if closeFn != nil {
			return closeFn()
		}
		return nil
	},
}

// gtSite is a requested variant position. ref/alt are empty unless supplied, in
// which case hasRA requires the matched record to share them.
type gtSite struct {
	chrom    string
	pos      int
	ref, alt string
	hasRA    bool
}

// gtRow is one output line: a matched record's genotype-class counts.
type gtRow struct {
	chrom  string
	pos    int
	ref    string
	alt    string
	counts map[string]int
}

// collectGtSites parses command-line loci (chrom:pos[:ref:alt]) followed by the
// sites in --sites FILE, preserving order.
func collectGtSites(loci []string, sitesFile string) ([]gtSite, error) {
	var sites []gtSite
	for _, locus := range loci {
		s, err := parseLocus(locus)
		if err != nil {
			return nil, err
		}
		sites = append(sites, s)
	}
	if sitesFile != "" {
		fileSites, err := readSitesFile(sitesFile)
		if err != nil {
			return nil, err
		}
		sites = append(sites, fileSites...)
	}
	return sites, nil
}

// parseLocus parses a "chrom:pos" or "chrom:pos:ref:alt" command-line locus.
func parseLocus(locus string) (gtSite, error) {
	parts := strings.Split(locus, ":")
	switch len(parts) {
	case 2:
		pos, err := strconv.Atoi(parts[1])
		if err != nil {
			return gtSite{}, fmt.Errorf("invalid locus %q: bad position", locus)
		}
		return gtSite{chrom: parts[0], pos: pos}, nil
	case 4:
		pos, err := strconv.Atoi(parts[1])
		if err != nil {
			return gtSite{}, fmt.Errorf("invalid locus %q: bad position", locus)
		}
		return gtSite{chrom: parts[0], pos: pos, ref: parts[2], alt: parts[3], hasRA: true}, nil
	default:
		return gtSite{}, fmt.Errorf("invalid locus %q: expected chrom:pos or chrom:pos:ref:alt", locus)
	}
}

// readSitesFile reads whitespace-separated "chrom pos [ref alt]" lines, skipping
// blank lines and '#' comments.
func readSitesFile(filename string) ([]gtSite, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var sites []gtSite
	scanner := bufio.NewScanner(f)
	line := 0
	for scanner.Scan() {
		line++
		text := strings.TrimSpace(scanner.Text())
		if text == "" || strings.HasPrefix(text, "#") {
			continue
		}
		fields := strings.Fields(text)
		pos, err := strconv.Atoi(fields[1])
		if len(fields) < 2 || err != nil {
			return nil, fmt.Errorf("%s:%d: expected 'chrom pos [ref alt]'", filename, line)
		}
		site := gtSite{chrom: fields[0], pos: pos}
		if len(fields) >= 4 {
			site.ref, site.alt, site.hasRA = fields[2], fields[3], true
		}
		sites = append(sites, site)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return sites, nil
}

// gtCountSite queries the indexed VCF for site and returns one gtRow per matched
// record (empty when nothing matches).
func gtCountSite(ir *vcf.IndexedVcfReader, site gtSite, passingOnly bool) ([]*gtRow, error) {
	if !ir.HasRef(site.chrom) {
		return nil, nil
	}
	seq, err := ir.Query(site.chrom, site.pos-1, site.pos)
	if err != nil {
		return nil, err
	}
	var rows []*gtRow
	for rec, qerr := range seq {
		if qerr != nil {
			return nil, qerr
		}
		if rec.Pos != site.pos {
			continue // a spanning record that merely overlaps the base
		}
		if passingOnly && rec.IsFiltered() {
			continue
		}
		if site.hasRA && (rec.Ref != site.ref || !containsStr(rec.Alt(), site.alt)) {
			continue
		}
		row, err := countGenotypes(rec)
		if err != nil {
			return nil, err
		}
		rows = append(rows, row)
	}
	return rows, nil
}

// countGenotypes tallies canonicalized GT classes across all samples of rec.
func countGenotypes(rec *vcf.VcfRecord) (*gtRow, error) {
	counts := map[string]int{}
	for i := 0; i < rec.NumSamples(); i++ {
		sample, err := rec.Sample(i)
		if err != nil {
			return nil, err
		}
		gt, ok := sample.Get("GT")
		if !ok {
			counts["."]++
			continue
		}
		counts[canonicalGT(gt.String())]++
	}
	return &gtRow{
		chrom:  rec.Chrom,
		pos:    rec.Pos,
		ref:    rec.Ref,
		alt:    orDot(rec.AltOrig()),
		counts: counts,
	}, nil
}

// canonicalGT collapses a raw GT string into an order-independent, unphased
// class: alleles are sorted (missing "." first) and joined with "/", so 1/0,
// 0/1 and 0|1 all become 0/1 and ./. stays ./.
func canonicalGT(raw string) string {
	if raw == "" || raw == "." {
		return "."
	}
	tokens := strings.FieldsFunc(raw, func(r rune) bool { return r == '/' || r == '|' })
	if len(tokens) == 0 {
		return "."
	}
	sort.Slice(tokens, func(i, j int) bool { return alleleRank(tokens[i]) < alleleRank(tokens[j]) })
	return strings.Join(tokens, "/")
}

// alleleRank orders allele tokens with missing ("." -> -1) first, then numeric,
// then any unrecognized token last.
func alleleRank(token string) int {
	if token == "." {
		return -1
	}
	n, err := strconv.Atoi(token)
	if err != nil {
		return math.MaxInt32
	}
	return n
}

// sortGtColumns orders genotype-class column names in place via gtColumnLess.
func sortGtColumns(cols []string) {
	sort.Slice(cols, func(i, j int) bool { return gtColumnLess(cols[i], cols[j]) })
}

// gtColumnLess orders genotype-class columns in the canonical VCF genotype order
// (grouped by the highest allele, then the next, ...): 0/0, 0/1, 1/1, 0/2, 1/2,
// 2/2, ... Missing-containing classes sort last, with ./. last of all. Tokens
// are already ascending (canonicalGT), so comparing from the last token down
// yields the highest-allele-first ordering.
func gtColumnLess(a, b string) bool {
	am, bm := strings.Contains(a, "."), strings.Contains(b, ".")
	if am != bm {
		return !am
	}
	at, bt := strings.Split(a, "/"), strings.Split(b, "/")
	for i := 1; i <= len(at) && i <= len(bt); i++ {
		ra, rb := columnRank(at[len(at)-i]), columnRank(bt[len(bt)-i])
		if ra != rb {
			return ra < rb
		}
	}
	return len(at) < len(bt)
}

// columnRank ranks allele tokens for column ordering, sorting missing (".")
// last (the opposite of alleleRank, which orders within a single genotype).
func columnRank(token string) int {
	if token == "." {
		return math.MaxInt32
	}
	n, err := strconv.Atoi(token)
	if err != nil {
		return math.MaxInt32 - 1
	}
	return n
}

func containsStr(values []string, target string) bool {
	for _, v := range values {
		if v == target {
			return true
		}
	}
	return false
}

func orDot(s string) string {
	if s == "" {
		return "."
	}
	return s
}

func init() {
	f := vcfGtCountCmd.Flags()
	f.StringVarP(&vcfGtCountOutput, "output", "o", "-", "Output filename (- for stdout)")
	f.StringVar(&vcfGtCountSites, "sites", "", "Read additional sites from FILE (chrom pos [ref alt] per line)")
	f.BoolVar(&vcfGtCountPassing, "passing", false, "Only count records that pass FILTER")
}
