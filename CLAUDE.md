# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CGLTK is a Go toolkit for computational genomics research. It provides CLI commands for sequence analysis, NGS data wrangling, and bioinformatics operations, with particular focus on Oxford Nanopore (long-read) sequence processing.

**Module:** `github.com/compgen-io/cgkit`
**Go version:** 1.23
**CLI framework:** Cobra

## Commands

```bash
# Build all targets (darwin_arm64, linux_arm64, linux_amd64)
make build

# Run all tests
make test
# equivalent to:
GOCACHE=/tmp/go-build-cache go test ./...

# Run a single test
go test ./align/... -run TestCigarCondense

# Run with CPU profiling
./cgkit --profile=cpu.prof <subcommand>
```

## Architecture

### Package Layout

- **`seqio/`** — FASTA/FASTQ readers with gzip support. Core type is `SeqQual`, which holds sequence, quality scores, name, and strand. Readers are streaming via `NextSeq()`.
- **`align/`** — Smith-Waterman local alignment with affine gap penalties. Includes special handling for Oxford Nanopore homopolymer error profiles.
- **`htsio/`** — SAM/BAM/CRAM reading and writing. Native BAM and SAM readers/writers; samtools only for CRAM. Includes BAI/TBI/CSI index parsers, tabix reader/writer, sorted BAM writer with merge sort.
- **`htsio/bgzf/`** — BGZF (Blocked GNU Zip Format) reader, writer, and indexed reader with LRU block cache. Used by BAM and tabix layers.
- **`support/sequtils/`** — DNA utilities: IUPAC ambiguity code matching, reverse complement, homopolymer run analysis, 4-bit DNA encoding.
- **`support/utils/`** — General utilities: semaphore for concurrency, float formatting, position-tracking reader.
- **`internal/cmd/`** — Cobra CLI commands organized into subpackages: `fastacmd`, `fastqcmd`, `seqcmd`, `ontcmd`, `samcmd`.

### HTS I/O System

The `htsio/` package provides native SAM/BAM I/O without external dependencies (samtools only for CRAM):

- `SamReader` interface: `Next()`, `Header()`, `Query()`, `Close()`
- `Query()` returns `iter.Seq2[*SamRecord, error]` — uses Go 1.23 range-over-func
- `NewSamReader()` auto-detects: `.bam` → `BamReader`, `.sam`/`.sam.gz` → `SamTextReader`, `.cram` → `SamtoolsSamReader`
- `NewSamWriter()` auto-selects: BAM (sorted/unsorted) → native, CRAM → samtools
- All query coordinates are 0-based half-open
- `ParseRegion()` converts samtools-style strings to 0-based half-open
- `IterReader()` bridges `iter.Seq2` back to `SamReader` for legacy callers
- `TabixReader`/`TabixWriter` handle tabix-indexed text files (BED, VCF, GFF) with TBI or CSI indexes
- `bgzf.IndexedReader` has an LRU block cache shared by BAI and tabix query paths

### Alignment System

The aligner (`align/`) is the most complex component:

- `NewLocalAligner()` — Smith-Waterman with soft clipping (for partial matches)
- `NewGlobalAligner()` — Full-sequence alignment
- `DnaAlignmentDefaults()` — Presets for Illumina short reads
- `OntAlignmentDefaults()` — Presets for Oxford Nanopore (looser gap penalties, homopolymer discounts)
- `AlignBatch()` — Parallel alignment using a semaphore-controlled goroutine pool
- Homopolymer discounts are precalculated and cached for performance

CIGAR strings use standard ops: M (match), I (insertion), D (deletion), S (soft clip). Helper functions `CigarCondense`/`CigarExpand` convert between run-length encoded and per-base forms.

### CLI Command Structure

Commands are registered in `internal/cmd/root.go` and grouped by file format or domain:
- `fasta-cat`, `fasta-wrap`, `fasta-gc` — FASTA operations
- `fastq-gc` — FASTQ operations
- `seq-pairwise`, `seq-revcomp` — Sequence analysis
- `ont-primers` — ONT primer detection/trimming with alignment statistics
- `ont-umi-dedup` — UMI deduplication: selects one representative per MI group from coordinate-sorted BAM. Secondary/supplementary alignments are dropped (cannot be reliably resolved in coordinate order). Supports `--threads` for parallel BGZF compression.
