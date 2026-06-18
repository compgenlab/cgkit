# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

cgio is a Go CLI toolkit for computational genomics research. It provides commands for sequence analysis, NGS data wrangling, and bioinformatics operations, with particular focus on Oxford Nanopore (long-read) sequence processing. The underlying library (sequence I/O, alignment, SAM/BAM/CRAM handling) lives in the separate `hts` module (`github.com/compgenlab/hts`).

**Module:** `github.com/compgenlab/cgio`
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
./cgio --profile=cpu.prof <subcommand>
```

## Dependency on hts

All format I/O and algorithms come from `github.com/compgenlab/hts` (packages
`seqio`, `align`, `htsio` and its subpackages, `support/*`, `analysis/seq`).
During local development the `hts` dependency is resolved through a `go.work`
workspace that joins a sibling `hts` checkout; release builds use the pinned
module version in `go.mod`. The `Makefile` deliberately does **not** set
`GOWORK=off`, so local builds pick up the workspace.

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
