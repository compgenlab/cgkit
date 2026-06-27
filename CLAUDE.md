# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

cgkit is a Go CLI toolkit for computational genomics research. It provides commands for sequence analysis, NGS data wrangling, and bioinformatics operations, with particular focus on Oxford Nanopore (long-read) sequence processing. The underlying library (sequence I/O, alignment, SAM/BAM/CRAM handling) lives in the separate `hts` module (`github.com/compgenlab/hts`).

**Module:** `github.com/compgenlab/cgkit`
**Go version:** 1.23
**CLI framework:** Cobra
**Library dependency:** `github.com/compgenlab/hts`

## Commands

```bash
# Build all targets (darwin_arm64, linux_arm64, linux_amd64)
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

## Dependency on hts

All format I/O and algorithms come from `github.com/compgenlab/hts` (packages
`seqio`, `align`, `htsio` and its subpackages, `support/*`, `analysis/seq`).

How the `hts` dependency resolves, by context:
- **Local builds** use the `go.work` workspace (parent dir, untracked) that joins
  a sibling `hts` checkout, so you build against your live local `hts` tree. The
  `Makefile` deliberately does **not** set `GOWORK=off`.
- **Remote/CI builds** (no `go.work` present) use the **latest released hts from
  GitHub**: the GitHub Actions workflow runs `go get github.com/compgenlab/hts@latest`
  before vet/test/build, with `GOPRIVATE=github.com/compgenlab/*` so a freshly
  pushed hts tag is fetched directly from GitHub (no module-proxy/sumdb lag).
- The committed `go.mod` pin is the fallback for `go install` users and the
  source archive; keep it current with `make bump-hts`.

### Cutting a release
The hts tag must land on GitHub before cgkit builds against it:
1. **hts**: tag `vX.Y.Z` on `main`, push the tag.
2. **cgkit**: `make bump-hts` (pins `go.mod` to the new hts), commit
   `go.mod`/`go.sum`, push. CI's `go get hts@latest` then resolves the same tag.

## Architecture

This repo holds only the CLI layer: `main.go` (entry point with `--profile`
support) and `internal/cmd/` (Cobra commands). The third-party dependencies are
cobra/pflag; everything genomics-related is delegated to `hts`.

### CLI Command Structure

Commands are registered in `internal/cmd/root.go` and grouped by file format or domain:
- `fasta-cat`, `fasta-wrap`, `fasta-gc` — FASTA operations
- `fastq-gc` — FASTQ operations
- `sam-stats` — Summary statistics for SAM/BAM/CRAM: read counts, mapping rates, Q30, depth, SAM flag breakdown, per-reference counts, optional `--tags` value distributions and `--calc-insert` median. Profiles the first read of each pair only (ports `ngsutils bam-stats`). Phase 1 omits the `--gtf` gene-model and `--bed` on-target stats.
- `seq-pairwise`, `seq-revcomp` — Sequence analysis
- `ont-primers` — ONT primer detection/trimming with alignment statistics
- `ont-umi-dedup` — UMI deduplication: selects one representative per MI group from coordinate-sorted BAM. Secondary/supplementary alignments are dropped (cannot be reliably resolved in coordinate order). Supports `--threads` for parallel BGZF compression.
