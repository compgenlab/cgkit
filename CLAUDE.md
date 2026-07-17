# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

cgkit is a Go CLI toolkit for computational genomics research. It provides commands for sequence analysis, NGS data wrangling, and bioinformatics operations, with particular focus on Oxford Nanopore (long-read) sequence processing. The underlying library (sequence I/O, alignment, SAM/BAM/CRAM handling) lives in the separate `cghts` module (`github.com/compgenlab/cghts`).

**Module:** `github.com/compgenlab/cgkit`
**Go version:** 1.23
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
cobra/pflag; everything genomics-related is delegated to `cghts`.

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
