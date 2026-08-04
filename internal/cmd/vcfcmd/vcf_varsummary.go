package vcfcmd

import (
	"bufio"
	"context"
	"encoding/json"
	"fmt"
	"strings"

	"github.com/compgenlab/cghts/varstore"
	"github.com/spf13/cobra"
)

var (
	vcfVarSummaryFormat  string
	vcfVarSummarySamples bool
	vcfVarSummarySites   bool
	vcfVarSummaryCounts  bool
	vcfVarSummaryOutput  string
	vcfVarSummaryVerbose bool
	vcfVarSummaryStore   string
)

var vcfVarSummaryCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": "v0.6.0"},
	Use:         "vcf-varsummary <input.vcf | store>",
	Short:       "Describe a genotype store or VCF: samples, contigs, provenance",
	Long: `Report what a variant source contains, without querying genotypes.

The input may be a VCF (plain or bgzipped) or a Parquet store written by
vcf-toparquet, local or remote. The backend is inferred from the path; override
with --store.

The default report never reads a single record. It is O(header) on both
backends -- instant against a 200 GB whole-genome VCF -- because a summary that
silently takes twenty minutes is worse than one that tells you what it will
cost. Everything requiring a pass over the data is opt-in:

  --samples    the sample roster, one per line               (free on both)
  --counts     record totals and the per-chromosome census   (store: free)
  --sites      every variant, as chrom pos ref alt ...       (store: streamed)

For a store, --counts is read straight out of the manifest and costs nothing.
For a VCF there is no such record and no index of "every variant", so both
--counts and --sites are a full pass over the file; that is reported under -v
rather than sprung on you.

What each backend can say differs, and the report says which it is rather than
inventing the missing half. A store carries a manifest: when it was written, by
what command, from which inputs, at what --min-dp, how many rows landed in each
member, and how many sites and calls each chromosome contributed. A VCF carries
none of that -- it is not a conversion, so there is nothing to record about one.
What a VCF does carry is its header: the sample roster, the ##contig
declarations, and whether a tabix index sits beside it.

The per-chromosome census is the useful part of a store's manifest. It counts
the rows that were written, not the inputs that were requested, so it is the one
field that can contradict the rest: a conversion that stopped after chr3 still
names all 22 inputs and declares all 22 contigs, because both are stamped before
the first record is read. A store missing its manifest cannot be queried at all,
and this command is how you find out why.

For an indexed VCF the contigs that actually carry records come from the tabix
index, which is the cheapest useful thing a VCF summary can report and needs no
scan. Exact record counts are not recoverable from an index, so they are not
estimated from it -- they are reported as requiring --counts.

  --samples          list the sample roster, one per line
  --sites            stream the variant catalog
  --counts           report record totals and the per-chromosome census
  --format F         text (default) or json; json emits a store's manifest
                     verbatim, so "| jq" works on it
  -o, --output FILE  write here instead of stdout
  -v, --verbose      note what is being read, and what a scan will cost`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if vcfVarSummaryFormat != "text" && vcfVarSummaryFormat != "json" {
			return fmt.Errorf("unknown --format %q (use text or json)", vcfVarSummaryFormat)
		}

		ctx := cmd.Context()
		store, err := openVarStore(ctx, args[0], vcfVarSummaryStore)
		if err != nil {
			return err
		}
		defer store.Close()

		w, closeFn, err := openOutput(cmd, vcfVarSummaryOutput)
		if err != nil {
			return err
		}
		out := bufio.NewWriter(w)

		if err := writeSummary(cmd, out, store, args[0]); err != nil {
			return err
		}
		if err := out.Flush(); err != nil {
			return err
		}
		if closeFn != nil {
			return closeFn()
		}
		return nil
	},
}

func writeSummary(cmd *cobra.Command, out *bufio.Writer, store varstore.Store, path string) error {
	ctx := cmd.Context()
	errOut := cmd.ErrOrStderr()

	if vcfVarSummaryFormat == "json" {
		return writeSummaryJSON(out, store, path)
	}

	samples, err := store.Samples()
	if err != nil {
		return err
	}

	if vcfVarSummarySamples {
		for _, s := range samples {
			fmt.Fprintln(out, s)
		}
		return nil
	}
	if vcfVarSummarySites {
		if _, isVcf := store.(*varstore.VcfStore); isVcf && vcfVarSummaryVerbose {
			fmt.Fprintln(errOut, "note: a VCF has no catalog to read, so --sites is a full pass over the file")
		}
		fmt.Fprintln(out, strings.Join([]string{
			"chrom", "pos", "ref", "alt", "ac", "an", "n_carriers", "n_called", "n_lowdp",
		}, "\t"))
		return store.Sites(func(s varstore.Site) bool {
			fmt.Fprintf(out, "%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n",
				s.Chrom, s.Pos, s.Ref, s.Alt, s.AC, s.AN, s.NCarriers, s.NCalled, s.NLowDP)
			return true
		})
	}

	switch s := store.(type) {
	case *varstore.ParquetStore:
		return summarizeStore(out, errOut, s, path, samples)
	case *varstore.VcfStore:
		return summarizeVcf(ctx, out, errOut, s, path, samples)
	}
	return fmt.Errorf("unknown backend %T", store)
}

func summarizeStore(out *bufio.Writer, errOut interface{ Write([]byte) (int, error) },
	s *varstore.ParquetStore, path string, samples []string) error {
	m := s.Manifest()
	p := s.Provenance()

	fmt.Fprintf(out, "source     parquet store %s\n", path)
	fmt.Fprintf(out, "samples    %d\n", len(samples))
	fmt.Fprintf(out, "created    %s\n", m.Created.Format("2006-01-02 15:04:05 MST"))
	if m.Program != "" {
		fmt.Fprintf(out, "written by %s\n", m.Program)
	}
	if len(m.Sources) > 0 {
		fmt.Fprintf(out, "converted from\n")
		for _, src := range m.Sources {
			fmt.Fprintf(out, "  %s\n", src)
		}
	}
	fmt.Fprintf(out, "spans      %s", p.Spans)
	if p.Spans == varstore.SpansSites {
		fmt.Fprintf(out, "  (answers only for variants in the sites catalog)")
	}
	fmt.Fprintln(out)
	if p.NoCallable {
		fmt.Fprintf(out, "callable   not tracked (--no-callable): --hom-ref will refuse\n")
	} else {
		fmt.Fprintf(out, "min-dp     %d at conversion\n", p.MinDP)
	}

	fmt.Fprintf(out, "\nmembers\n")
	for _, name := range []string{
		varstore.CallsMember, varstore.SitesMember, varstore.RegionsMember,
	} {
		info, ok := m.Members[name]
		if !ok {
			continue
		}
		fmt.Fprintf(out, "  %-9s %12d rows  %14d bytes\n", name, info.Rows, info.Bytes)
	}

	if !vcfVarSummaryCounts {
		return nil
	}
	fmt.Fprintf(out, "\nper chromosome (as written, not as requested)\n")
	fmt.Fprintf(out, "  %-12s %12s %12s  %s\n", "chrom", "sites", "calls", "range")
	for _, c := range m.Chromosomes {
		fmt.Fprintf(out, "  %-12s %12d %12d  %d-%d\n", c.Name, c.Sites, c.Calls, c.FirstPos, c.LastPos)
	}
	// A contig declared but never written is what a stopped conversion looks
	// like, and it is invisible without saying so: the declaration is stamped
	// before the first record is read.
	written := map[string]bool{}
	for _, c := range m.Chromosomes {
		written[c.Name] = true
	}
	var missing []string
	for _, id := range contigIDs(m.ContigsDeclared) {
		if !written[id] {
			missing = append(missing, id)
		}
	}
	if len(missing) > 0 {
		fmt.Fprintf(out, "\n  NOTE: declared by the source but holding no sites: %s\n",
			strings.Join(missing, ", "))
		fmt.Fprintf(out, "        that is normal for a subset conversion, and is what an\n")
		fmt.Fprintf(out, "        interrupted one would also look like.\n")
	}
	return nil
}

func summarizeVcf(ctx context.Context, out *bufio.Writer, errOut interface{ Write([]byte) (int, error) },
	s *varstore.VcfStore, path string, samples []string) error {
	fmt.Fprintf(out, "source     vcf %s\n", path)
	fmt.Fprintf(out, "samples    %d\n", len(samples))
	if s.Indexed() {
		fmt.Fprintf(out, "indexed    yes\n")
	} else {
		fmt.Fprintf(out, "indexed    no  (region queries scan the whole file)\n")
	}
	fmt.Fprintf(out, "provenance none -- a VCF is not a conversion, so there is nothing recorded\n")
	fmt.Fprintf(out, "           about how it was made; --min-dp and span semantics are\n")
	fmt.Fprintf(out, "           properties of a store, not of this file\n")

	if !vcfVarSummaryCounts {
		fmt.Fprintf(out, "\ncounts     require --counts (a full pass over the file)\n")
		return nil
	}

	if vcfVarSummaryVerbose {
		fmt.Fprintln(errOut, "note: counting requires reading every record")
	}
	type chromCount struct {
		name              string
		sites             int64
		firstPos, lastPos int32
	}
	var order []*chromCount
	idx := map[string]*chromCount{}
	var total int64
	if err := s.Sites(func(site varstore.Site) bool {
		total++
		c, ok := idx[site.Chrom]
		if !ok {
			c = &chromCount{name: site.Chrom, firstPos: site.Pos, lastPos: site.Pos}
			idx[site.Chrom] = c
			order = append(order, c)
		}
		c.sites++
		if site.Pos < c.firstPos {
			c.firstPos = site.Pos
		}
		if site.Pos > c.lastPos {
			c.lastPos = site.Pos
		}
		return true
	}); err != nil {
		return err
	}

	fmt.Fprintf(out, "\ncounts\n")
	fmt.Fprintf(out, "  sites    %d  (one per ALT allele)\n", total)
	fmt.Fprintf(out, "\nper chromosome\n")
	fmt.Fprintf(out, "  %-12s %12s  %s\n", "chrom", "sites", "range")
	for _, c := range order {
		fmt.Fprintf(out, "  %-12s %12d  %d-%d\n", c.name, c.sites, c.firstPos, c.lastPos)
	}
	return nil
}

// writeSummaryJSON emits a store's manifest verbatim, so that gzipping it costs
// nothing in scriptability: "| jq" is one pipe away rather than zero.
func writeSummaryJSON(out *bufio.Writer, store varstore.Store, path string) error {
	enc := json.NewEncoder(out)
	enc.SetIndent("", "  ")
	switch s := store.(type) {
	case *varstore.ParquetStore:
		return enc.Encode(s.Manifest())
	case *varstore.VcfStore:
		samples, err := s.Samples()
		if err != nil {
			return err
		}
		return enc.Encode(struct {
			Source  string   `json:"source"`
			Kind    string   `json:"kind"`
			Indexed bool     `json:"indexed"`
			Samples []string `json:"samples"`
		}{path, "vcf", s.Indexed(), samples})
	}
	return fmt.Errorf("unknown backend %T", store)
}

// contigIDs pulls the ID out of each ##contig line.
func contigIDs(lines []string) []string {
	var out []string
	for _, l := range lines {
		i := strings.Index(l, "ID=")
		if i < 0 {
			continue
		}
		rest := l[i+3:]
		if j := strings.IndexAny(rest, ",>"); j >= 0 {
			rest = rest[:j]
		}
		if rest != "" {
			out = append(out, rest)
		}
	}
	return out
}

func init() {
	f := vcfVarSummaryCmd.Flags()
	f.BoolVar(&vcfVarSummarySamples, "samples", false, "List the sample roster, one per line")
	f.BoolVar(&vcfVarSummarySites, "sites", false, "Stream the variant catalog (a full pass for a VCF)")
	f.BoolVar(&vcfVarSummaryCounts, "counts", false, "Report totals and the per-chromosome census (a full pass for a VCF)")
	f.StringVar(&vcfVarSummaryFormat, "format", "text", "Output format: text or json")
	f.StringVarP(&vcfVarSummaryOutput, "output", "o", "-", "Output filename")
	f.BoolVarP(&vcfVarSummaryVerbose, "verbose", "v", false, "Note what is being read, on stderr")
	f.StringVar(&vcfVarSummaryStore, "store", "", "Force the backend: vcf or parquet")
}
