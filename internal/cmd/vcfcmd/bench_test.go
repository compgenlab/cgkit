package vcfcmd

import (
	"bufio"
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/varstore"
)

// Benchmarks for genotype queries, over a corpus this file GENERATES rather than a
// fixture it commits.
//
// That is deliberate. Every performance number this project has quoted so far
// lives in a commit message and was measured on 1000 Genomes data that is not in
// the repository, so none of it can be reproduced or re-checked. A seeded
// generator can: `go test -bench . ./internal/cmd/vcfcmd/` rebuilds the corpus and
// re-measures.
//
// Two questions are asked. First, whether one bulk query beats looping a
// single-locus query -- the reason the Store API collapsed to one Calls method,
// and so far an arithmetic argument rather than a measurement. Second, how the
// Parquet store compares with querying the plain VCF it came from, which
// docs/vcf-tovarstore.md currently states outright has never been measured.
//
// Wall time and row counts only. Bytes read and row groups decoded would be the
// more durable numbers, but they need a counter inside varstore, which lives in
// cghts.

// runVcfTB is runVcf for benchmarks, which get a testing.B rather than a T.
func runVcfTB(tb testing.TB, args ...string) string {
	tb.Helper()
	root, buf := vcfTestRoot(args...)
	if err := root.Execute(); err != nil {
		tb.Fatalf("Execute(%v): %v", args, err)
	}
	return buf.String()
}

// genCallset writes a synthetic VCF and converts it to a store, returning both.
//
// Genotypes are mostly reference, as a real cohort callset is -- that sparsity is
// the whole premise of the store format, so a corpus without it would flatter the
// store by making both files the same size.
func genCallset(tb testing.TB, nSamples, nSites, rowGroup int, seed int64) (vcfPath, storeBase string) {
	tb.Helper()
	dir := tb.TempDir()
	vcfPath = filepath.Join(dir, "synth.vcf")

	f, err := os.Create(vcfPath)
	if err != nil {
		tb.Fatal(err)
	}
	w := bufio.NewWriterSize(f, 1<<20)

	samples := make([]string, nSamples)
	for i := range samples {
		samples[i] = fmt.Sprintf("S%04d", i)
	}
	fmt.Fprintln(w, "##fileformat=VCFv4.2")
	fmt.Fprintf(w, "##contig=<ID=chr1,length=%d>\n", nSites*10+1000)
	fmt.Fprintln(w, `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`)
	fmt.Fprintln(w, `##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">`)
	fmt.Fprintln(w, `##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">`)
	fmt.Fprintln(w, "#"+strings.Join(append([]string{
		"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
	}, samples...), "\t"))

	rng := rand.New(rand.NewSource(seed))
	bases := []string{"A", "C", "G", "T"}
	var row strings.Builder
	for s := 0; s < nSites; s++ {
		pos := s*10 + 100
		ref := bases[rng.Intn(4)]
		alt := bases[(rng.Intn(3)+1+indexOf(bases, ref))%4]
		row.Reset()
		fmt.Fprintf(&row, "chr1\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT:DP:GQ", pos, ref, alt)
		for i := 0; i < nSamples; i++ {
			// ~5% of genotypes carry the alternate, so ~95% of the matrix is
			// reference and the store stays sparse.
			gt := "0/0"
			switch n := rng.Float64(); {
			case n < 0.045:
				gt = "0/1"
			case n < 0.05:
				gt = "1/1"
			}
			fmt.Fprintf(&row, "\t%s:%d:99", gt, 20+rng.Intn(30))
		}
		fmt.Fprintln(w, row.String())
	}
	if err := w.Flush(); err != nil {
		tb.Fatal(err)
	}
	if err := f.Close(); err != nil {
		tb.Fatal(err)
	}

	storeBase = filepath.Join(dir, "store") + string(os.PathSeparator)
	runVcfTB(tb, "vcf-tovarstore", "--out", storeBase,
		"--row-group-size", strconv.Itoa(rowGroup), vcfPath)
	return vcfPath, storeBase
}

func indexOf(ss []string, s string) int {
	for i, v := range ss {
		if v == s {
			return i
		}
	}
	return 0
}

// panelOf picks nTargets loci spread across the callset, as a variant panel is.
func panelOf(tb testing.TB, storeBase string, nTargets int) []varstore.Locus {
	tb.Helper()
	s, err := varstore.OpenParquet(storeBase)
	if err != nil {
		tb.Fatal(err)
	}
	defer s.Close()

	var all []varstore.Locus
	if err := s.Sites(func(site varstore.Site) bool {
		all = append(all, site.Locus())
		return true
	}); err != nil {
		tb.Fatal(err)
	}
	if nTargets > len(all) {
		nTargets = len(all)
	}
	step := len(all) / nTargets
	out := make([]varstore.Locus, 0, nTargets)
	for i := 0; i < nTargets; i++ {
		out = append(out, all[i*step])
	}
	return out
}

// benchSizes returns the corpus dimensions, overridable so a run can be scaled
// without editing this file.
//
// rowGroup matters more than it looks. Left at the converter's 250,000 default a
// modest corpus lands in ONE row group, and then position statistics can never
// exclude it -- every locus lookup scans the whole file, which is the loop's worst
// case rather than a realistic one. Real stores hold many groups, so the default
// here is small enough to produce several.
func benchSizes() (nSamples, nSites, rowGroup int) {
	nSamples, nSites, rowGroup = 200, 20000, 10000
	if os.Getenv("CGKIT_BENCH_BIG") != "" {
		nSamples, nSites = 500, 50000
	}
	atoiEnv := func(key string, into *int) {
		if v := os.Getenv(key); v != "" {
			if n, err := strconv.Atoi(v); err == nil && n > 0 {
				*into = n
			}
		}
	}
	atoiEnv("CGKIT_BENCH_SAMPLES", &nSamples)
	atoiEnv("CGKIT_BENCH_SITES", &nSites)
	atoiEnv("CGKIT_BENCH_ROWGROUP", &rowGroup)
	return nSamples, nSites, rowGroup
}

// BenchmarkQueryPanel is the measurement the single-Calls API was argued for: one
// bulk query against one query per locus, over the same panel.
func BenchmarkQueryPanel(b *testing.B) {
	nSamples, nSites, rowGroup := benchSizes()
	_, storeBase := genCallset(b, nSamples, nSites, rowGroup, 1)

	store, err := varstore.OpenParquet(storeBase)
	if err != nil {
		b.Fatal(err)
	}
	defer store.Close()

	for _, nTargets := range []int{1, 10, 100, 1000} {
		panel := panelOf(b, storeBase, nTargets)

		b.Run(fmt.Sprintf("bulk/%d", nTargets), func(b *testing.B) {
			for b.Loop() {
				rows, err := varstore.CollectCalls(store, varstore.Query{Loci: panel})
				if err != nil {
					b.Fatal(err)
				}
				b.ReportMetric(float64(len(rows)), "rows")
			}
		})

		b.Run(fmt.Sprintf("loop/%d", nTargets), func(b *testing.B) {
			for b.Loop() {
				n := 0
				for _, l := range panel {
					rows, err := varstore.CollectCalls(store,
						varstore.Query{Loci: []varstore.Locus{l}})
					if err != nil {
						b.Fatal(err)
					}
					n += len(rows)
				}
				b.ReportMetric(float64(n), "rows")
			}
		})
	}
}

// BenchmarkStoreVsVcf is the comparison docs/vcf-tovarstore.md says has never been
// made: the same panel query against the store and against the VCF it came from.
func BenchmarkStoreVsVcf(b *testing.B) {
	nSamples, nSites, rowGroup := benchSizes()
	vcfPath, storeBase := genCallset(b, nSamples, nSites, rowGroup, 1)

	store, err := varstore.OpenParquet(storeBase)
	if err != nil {
		b.Fatal(err)
	}
	defer store.Close()
	plain, err := varstore.OpenVcf(vcfPath)
	if err != nil {
		b.Fatal(err)
	}
	defer plain.Close()

	// Also measure a bgzipped, tabix-indexed VCF. Without this the comparison
	// flatters the store: an unindexed VCF has to be scanned end to end for any
	// lookup, which is a property of the file not being indexed rather than of the
	// format. With an index the VCF backend can seek.
	gzPath := vcfPath + ".gz"
	runVcfTB(b, "vcf-filter", "-o", gzPath, "--tbi", vcfPath)
	indexed, err := varstore.OpenVcf(gzPath)
	if err != nil {
		b.Fatal(err)
	}
	defer indexed.Close()

	if st, err := os.Stat(vcfPath); err == nil {
		b.Logf("corpus: %d samples x %d sites, VCF %d bytes", nSamples, nSites, st.Size())
	}

	for _, nTargets := range []int{1, 100} {
		panel := panelOf(b, storeBase, nTargets)
		q := varstore.Query{Loci: panel}
		for _, tc := range []struct {
			name string
			s    varstore.Store
		}{{"store", store}, {"vcf", plain}, {"vcf-indexed", indexed}} {
			b.Run(fmt.Sprintf("%s/%d", tc.name, nTargets), func(b *testing.B) {
				for b.Loop() {
					if _, err := varstore.CollectCalls(tc.s, q); err != nil {
						b.Fatal(err)
					}
				}
			})
		}
	}
}
