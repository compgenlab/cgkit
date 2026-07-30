# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

cgkit is a Go CLI toolkit for computational genomics research. It provides commands for sequence analysis, NGS data wrangling, and bioinformatics operations, with particular focus on Oxford Nanopore (long-read) sequence processing. The underlying library (sequence I/O, alignment, SAM/BAM/CRAM handling) lives in the separate `cghts` module (`github.com/compgenlab/cghts`).

**Module:** `github.com/compgenlab/cgkit`
**Go version:** 1.24.9
**CLI framework:** Cobra
**Library dependency:** `github.com/compgenlab/cghts`

## Commands

```bash
# Build all targets (darwin_arm64, darwin_amd64, linux_arm64, linux_amd64, windows_amd64)
make build

# Run all tests
make test
# equivalent to:
GOCACHE=/tmp/go-build-cache go test ./...

# Run a single test
go test ./internal/cmd/samcmd/... -run TestSamStats

# Run with CPU profiling
./cgkit --profile=cpu.prof <subcommand>
```

## Dependency on cghts

All format I/O and algorithms come from `github.com/compgenlab/cghts` (packages
`seqio`, `align`, `htsio` and its subpackages, `support/*`, `analysis/seq`).

How the `cghts` dependency resolves, by context:
- **Local builds** use the `go.work` workspace (parent dir, untracked) that joins
  a sibling `cghts` checkout, so you build against your live local `cghts` tree. The
  `Makefile` deliberately does **not** set `GOWORK=off`.
- **Remote/CI builds** (no `go.work` present) use the **latest released cghts from
  GitHub**: the GitHub Actions workflow runs `go get github.com/compgenlab/cghts@latest`
  before vet/test/build, with `GOPRIVATE=github.com/compgenlab/*` so a freshly
  pushed cghts tag is fetched directly from GitHub (no module-proxy/sumdb lag).
- The committed `go.mod` pin is the fallback for `go install` users and the
  source archive; keep it current with `make bump-cghts`.

### Cutting a release
The cghts tag must land on GitHub before cgkit builds against it:
1. **cghts**: tag `vX.Y.Z` on `main`, push the tag.
2. **cgkit**: `make bump-cghts` (pins `go.mod` to the new cghts), commit
   `go.mod`/`go.sum`, push. CI's `go get cghts@latest` then resolves the same tag.

## Architecture

This repo holds only the CLI layer: `main.go` (entry point with `--profile`
support) and `internal/cmd/` (Cobra commands). The third-party dependencies are
cobra/pflag plus `parquet-go/parquet-go`; everything genomics-related is
delegated to `cghts`.

`internal/varstore/` is the one exception to "CLI layer only": it holds the
on-disk schema for the Parquet genotype store and a `Store` interface with VCF
and Parquet implementations. It lives here rather than in `cghts` because it is
storage/IO glue rather than a genomics algorithm. Its `vcfrecord.go` is the
single authoritative reading of a VCF genotype — both `vcf-toparquet` and the
VCF-backed `Store` go through it, so a query against a VCF and the same query
against a store converted from it cannot drift apart.

### CLI Command Structure

Commands are registered in `internal/cmd/root.go` and grouped by file format or domain:
- `fasta-cat`, `fasta-wrap`, `fasta-gc` — FASTA operations
- `fastq-gc` — FASTQ operations
- `sam-stats` — Summary statistics for SAM/BAM/CRAM: read counts, mapping rates, Q30, depth, SAM flag breakdown, per-reference counts, optional `--tags` value distributions and `--calc-insert` median. Profiles the first read of each pair only (ports `ngsutils bam-stats`). Phase 1 omits the `--gtf` gene-model and `--bed` on-target stats.
- `seq-pairwise`, `seq-revcomp` — Sequence analysis
- `ont-polya` — Per-read poly(A)/cleavage site calling from a strand-specific aligned BAM. Finds the mRNA 3' end (FLAG 0x10, or `--antisense`), traces back through the tail with a windowed A-fraction test, and reports the first tail base's 1-based genomic position. The trace deliberately continues past the soft-clip boundary into aligned bases, since aligners absorb genome-encoded A's at real sites — which also makes the tool prone to reporting internal priming; `polya_source` (`--polya-src`) is the partial hook for filtering that downstream. Secondary/supplementary alignments are skipped. PAS motif annotation is deliberately out of scope: it is a per-site property needing a reference, so it belongs after clustering reads into sites.
- `ont-tags` — ONT adapter/primer detection and trimming from FASTQ, with alignment statistics (embeds a default primer set; override with `--primers-fasta`)
- `ont-umi-cluster` — Collapse similar UMIs in a coordinate-sorted BAM into `MI` groups
- `ont-umi-dedup` — UMI deduplication: selects one representative per MI group from coordinate-sorted BAM. Secondary/supplementary alignments are dropped (cannot be reliably resolved in coordinate order). Supports `--threads` for parallel BGZF compression.
- `ont-umi-lookup` — Match reads in an aligned BAM to UMI clusters from `ont-umi-cluster` output
- `vcf-toparquet` — Convert one or more VCFs into a sparse Parquet genotype store (several inputs is how whole-genome callsets ship: one per chromosome). Inputs must carry **exactly the same samples**; differing column order is **remapped**, because genotype columns are positional and getting that wrong silently attributes every genotype to the wrong person. A set mismatch errors, naming what differs. Inputs must not overlap (no revisiting a chromosome, no going backwards within one) since that would duplicate sites and split AC/AN; supplying them in coordinate order is not required for correctness but keeps pruning tight (~1.8x on a locus query). A failed conversion calls `Writer.Discard()` so no partial store is left to look valid or block the retry through the overwrite guard.: `BASE.calls.parquet` (one row per ALT-carrying genotype), `BASE.sites.parquet` (one row per interrogated site, carrying **AC/AN** allele counts alongside `n_carriers`/`n_called`/`n_lowdp` sample counts — these are not interchangeable: a 1/1 genotype is one carrier but two alt alleles, so AC ≥ n_carriers whenever a homozygote is present, and AN counts alleles ungated while n_called counts samples passing `--min-dp`. AC/AN are recomputed from the store's own samples rather than copied from source INFO, which would be wrong after splitting multiallelics or converting a cohort subset, and they are computed outside the `--no-callable` guard since allele counts are a property of genotypes, not coverage), `BASE.regions.parquet` (contiguous runs of covered sites per sample). All three are one inseparable set — and `--out` ending in `/` names a **directory** instead of a filename prefix, creating it (with parents) and writing bare `calls.parquet`/`sites.parquet`/`regions.parquet` inside, so the set is one thing to move. `vcf-varquery` accepts either form, plus the bare directory without the slash and any single member path; `varstore.MemberPath`/`TrimStoreSuffix` own that resolution. Conversion **refuses to overwrite an existing store** — any one of the three members already present under `--out`, or a prefix-form base naming an existing directory, stops it and asks for `--force`; writing truncates all three, and a half-replaced set is worse than either whole outcome. The guard keys on the members, so an existing directory holding unrelated files is fine and its contents are left alone. The sites file is **not** derivable from the calls: taking distinct loci out of the calls only recovers every site when the store holds an entire joint callset, and over a sample subset the sites nobody in that subset carries vanish silently, so a later query reports "never interrogated" for a position that was in fact observed and reference. Records are normalized to one variant per row (multiallelics split, focal allele recoded to 1, other alternates masked to `.` so a `1/2` sample is correctly a carrier of both; AD taken per allele, not summed). Indels are **not** left-aligned. A site counts as called only when the caller actually made a call there *and* DP ≥ `--min-dp` — `./.` at high depth is a declined call, not a covered one. **The regions file records runs of *catalog sites* at which a sample was called; the interval form is a compression of that per-site fact and claims nothing about the bases in between.** A plain VCF reports variants and asserts nothing about any other position, so the sites catalog is the exact boundary of what is answerable: an off-catalog locus is `not_assayed` for every sample, never a set of reference calls, even where run intervals bracket it (`Classify` returns early on the catalog check for precisely this reason, and `vcf-varquery` warns). Stores record `cgkit.span_semantics` = `sites`; only a gVCF-derived store (`blocks`, not yet produced) could answer off-catalog. See the GVCF memory. `-v/--verbose` reports progress, a per-chromosome trace, and a conversion summary on stderr — notably a census of which of DP/GQ/AD the input actually carries, since a gate can only act on a field the data has.
- `vcf-varquery` — Query which subjects carry a variant (`--variant chrom:pos:ref:alt`) or which variants a subject carries (`--sample`), over either a VCF or a `vcf-toparquet` store, chosen by path or `--store`. **Chromosome naming is auto-converted**: UCSC (`chr22`), Ensembl (`22`) and NCBI RefSeq (`NC_000022.11`) spellings all resolve, via cghts's `CanonicalContig`/`ContigConverter`. Comparisons use canonical identity, and tabix lookups resolve against the index's own `RefNames()` — tabix matches names verbatim, so a bare `22` against a `chr22` index used to fail with `unknown reference`. Conversion preserves the source's own spelling in the store rather than rewriting it. A contig the file lacks is an absence (no rows plus a warning), independent of whether the file is indexed; an unresolvable `--region` is an error, since that names a contig the caller asserts exists. `-v/--verbose` reports the chosen backend, the store's conversion provenance, per-variant AC/AN/AF and state tallies, and how many calls the quality gate actually excluded — including a warning when the gate could not act because the field is absent, and when `--min-dp` differs from the value the store's runs were built at. All verbose output goes to stderr so the tabular stream stays parseable. **Queries prune row groups rather than scanning the whole file** (`internal/varstore/prune.go`): input is coordinate-sorted and written in stream order, so parquet's per-group min/max on `pos` are tight, and a locus lookup skips every group that cannot contain it — measured 3.05x on a two-chromosome 1KG store (479ms to 157ms). Filters are conservative in one direction: they may keep a group holding nothing, never skip one holding a match, and the per-row checks are unchanged, so pruning can only affect how much is read. Chromosome bounds are only used when a group holds exactly one chromosome, since the statistics are lexicographic while chromosome names are not ordered that way. Note `--sample` queries cannot be pruned this way: every sample appears in nearly every row group, so the sample_id bloom filter always says maybe — making those fast would need a sample-major sort, which conflicts with the (chrom,pos) sort. `--classify` resolves every sample to carrier / uncertain / non-carrier / not-assayed rather than listing only carriers. The backends must agree: verified identical on 3,202 1KG samples, including at the worst-covered site in BRCA1. For that equivalence the query `--min-dp` should match the conversion `--min-dp`, since the Parquet side baked that threshold into its callable regions. An incomplete store (missing sites/regions, or built with `--no-callable`) makes `--classify` return `ErrNotClassifiable` rather than degrade — reporting an unobserved sample as a non-carrier is exactly the error the four states exist to prevent. `--hom-ref` adds **reference calls** to either mode's normal output rather than replacing it (the `gt` column separates them): in `--variant` mode the samples that are 0/0 at the locus, in `--sample` mode the catalog sites that subject was 0/0 at. It is **stricter than `--classify`'s `non_carrier`**, and the difference is the whole point: at a multiallelic record a 0/2 sample is not a carrier of allele 1, but it is not reference either, so it appears under neither — emitting `0/0` for it would be a genotype the source never contained. Both backends therefore key the exclusion on the source *record* (`Locus.Record()`, i.e. chrom+pos+ref), not the split locus, and the disqualifying ALT call is looked up ungated, since a below-gate ALT makes a sample an uncertain carrier rather than a reference observation. Everything that bounds `--classify` bounds this too: the gate must admit the reference call (a 0/0 at DP 3 under `--min-dp 10` is not an observation worth making), an off-catalog locus yields nothing, and an incomplete store returns `ErrNotClassifiable` in both modes. One asymmetry is unavoidable and is the one place the backends do not agree byte-for-byte: a store keeps only ALT genotypes, so a reference call recovered from one is a synthesized `0/0` with no DP/AD/GQ and no recoverable ploidy or phasing, where the same query against a VCF reports the recorded genotype and its quality fields — *which* rows appear is identical, and `-v` says so. `--format list` is refused with `--hom-ref` (bare sample ids cannot say which are carriers), as is `--classify` (it already resolves every sample). In `--sample` mode this walks the whole sites catalog — roughly one row per variant in the source — so `--region` is what keeps it bounded. Note `--format vcf` sorts its rows: loci are gathered per sample, so without that the second sample's private loci all land after the first's, and an unsorted VCF cannot be indexed.
