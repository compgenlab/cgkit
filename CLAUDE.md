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
- `vcf-toparquet` — Convert a VCF into a sparse Parquet genotype store: `BASE.calls.parquet` (one row per ALT-carrying genotype), `BASE.sites.parquet` (one row per interrogated site), `BASE.regions.parquet` (contiguous runs of covered sites per sample). All three are one inseparable set. The sites file is **not** derivable from the calls: taking distinct loci out of the calls only recovers every site when the store holds an entire joint callset, and over a sample subset the sites nobody in that subset carries vanish silently, so a later query reports "never interrogated" for a position that was in fact observed and reference. Records are normalized to one variant per row (multiallelics split, focal allele recoded to 1, other alternates masked to `.` so a `1/2` sample is correctly a carrier of both; AD taken per allele, not summed). Indels are **not** left-aligned. A site counts as called only when the caller actually made a call there *and* DP ≥ `--min-dp` — `./.` at high depth is a declined call, not a covered one. **The regions file records runs of *catalog sites* at which a sample was called; the interval form is a compression of that per-site fact and claims nothing about the bases in between.** A plain VCF reports variants and asserts nothing about any other position, so the sites catalog is the exact boundary of what is answerable: an off-catalog locus is `not_assayed` for every sample, never a set of reference calls, even where run intervals bracket it (`Classify` returns early on the catalog check for precisely this reason, and `vcf-varquery` warns). Stores record `cgkit.span_semantics` = `sites`; only a gVCF-derived store (`blocks`, not yet produced) could answer off-catalog. See the GVCF memory.
- `vcf-varquery` — Query which subjects carry a variant (`--variant chrom:pos:ref:alt`) or which variants a subject carries (`--sample`), over either a VCF or a `vcf-toparquet` store, chosen by path or `--store`. **Chromosome naming is auto-converted**: UCSC (`chr22`), Ensembl (`22`) and NCBI RefSeq (`NC_000022.11`) spellings all resolve, via cghts's `CanonicalContig`/`ContigConverter`. Comparisons use canonical identity, and tabix lookups resolve against the index's own `RefNames()` — tabix matches names verbatim, so a bare `22` against a `chr22` index used to fail with `unknown reference`. Conversion preserves the source's own spelling in the store rather than rewriting it. A contig the file lacks is an absence (no rows plus a warning), independent of whether the file is indexed; an unresolvable `--region` is an error, since that names a contig the caller asserts exists. `--classify` resolves every sample to carrier / uncertain / non-carrier / not-assayed rather than listing only carriers. The backends must agree: verified identical on 3,202 1KG samples, including at the worst-covered site in BRCA1. For that equivalence the query `--min-dp` should match the conversion `--min-dp`, since the Parquet side baked that threshold into its callable regions. An incomplete store (missing sites/regions, or built with `--no-callable`) makes `--classify` return `ErrNotClassifiable` rather than degrade — reporting an unobserved sample as a non-carrier is exactly the error the four states exist to prevent.
