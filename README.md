# cgltk

A Go toolkit for computational genomics research. Provides library packages for sequence I/O, alignment, and NGS data processing, along with CLI commands for common bioinformatics operations. Particular focus on Oxford Nanopore (long-read) sequencing workflows.

**Module:** `github.com/compgen-io/cgltk`

## Building

```bash
make build     # Build all targets (darwin_arm64, linux_arm64, linux_amd64)
make test      # Run all tests
```

## Library packages

### seqio — FASTA/FASTQ I/O

Streaming readers and writers for FASTA and FASTQ files with transparent gzip support.

- `SeqReader` / `SeqRecord` interfaces for uniform access across formats
- `FastaReader` / `FastqReader` — lazy, streaming readers via `NextSeq()`; support indexed lookup by name
- `FastaWriter` / `FastqWriter` — writers with optional line wrapping (FASTA) and gzip output
- `SeqQual` — core type holding sequence, quality, name, strand, and position; supports `RevComp()` and `Sub()` extraction
- Memory-efficient chunked iteration via Go `iter.Seq`

### align — Pairwise and multiple sequence alignment

Smith-Waterman based alignment with affine gap penalties and Oxford Nanopore-aware homopolymer discounting.

- `NewLocalAligner()` — Smith-Waterman local alignment (soft clipping)
- `NewGlobalAligner()` — Needleman-Wunsch end-to-end alignment
- `NewSemiGlobalAligner()` — full query aligned, free target end gaps
- `DnaAlignmentDefaults()` / `OntAlignmentDefaults()` — preset scoring parameters
- Configurable scoring matrix, gap penalties, clipping, and homopolymer discount via builder pattern
- `AlignBatch()` — parallel alignment with semaphore-controlled goroutine pool
- `CigarCondense()` / `CigarExpand()` — convert between run-length and per-base CIGAR formats
- `MSA()` — incremental consensus multiple sequence alignment with `Profile` output

### htsio — SAM/BAM/CRAM I/O

Reading and writing alignment files via samtools integration.

- `SamReader` — read SAM/BAM/CRAM with region, flag, MAPQ, and tag filters
- `SamWriter` — write SAM, BAM, or CRAM with optional coordinate/name sorting and multithreaded compression
- `SamRecord` — full SAM record with flag accessors (`IsUnmapped()`, `IsReverse()`, etc.) and typed tag access
- `SamHeader` — header manipulation including `@PG` line generation with auto-versioning
- `TagFilter` — flexible tag-based filtering with comparison operators

### tabix — BGZF compression

- `BGZipWriter` — block-gzipped output compatible with tabix indexing

### support packages

- **sequtils** — IUPAC ambiguity matching, reverse complement, homopolymer run analysis, 4-bit DNA encoding
- **utils** — `Semaphore` for concurrency control, `PositionTrackingReader`, float formatting
- **analysis/seq** — GC content calculation

## CLI commands

Usage: `cgltk [--profile=cpu.prof] <command>`

### FASTA

| Command | Description |
|---------|-------------|
| `fasta-wrap` | Reformat sequences to a specified line width (`-w`, default 70) |
| `fasta-gc` | Calculate GC content of sequences |

### FASTQ

| Command | Description |
|---------|-------------|
| `fastq-gc` | Calculate GC content of sequences |
| `fastq-tag` | Add a tag to the comment field of records |

### Sequence

| Command | Description |
|---------|-------------|
| `seq-revcomp` | Reverse complement a sequence |
| `seq-pairwise` | Pairwise alignment with configurable scoring, gap penalties, and homopolymer discounts |
| `seq-consensus-msa` | Multiple sequence alignment via incremental consensus (`--ont` for Nanopore defaults) |

### SAM/BAM/CRAM

| Command | Description |
|---------|-------------|
| `sam-export` | Export selected columns and tags as tab-delimited text |
| `sam-filter` | Filter reads by region, flags, MAPQ, or tags and write to a new file |
| `sam-toseq` | Convert reads to FASTA or FASTQ |

### Oxford Nanopore

| Command | Description |
|---------|-------------|
| `ont-tags` | Find and trim common ONT adapter/primer tags from FASTQ reads |
| `ont-umi-cluster` | Collapse similar UMIs in a coordinate-sorted BAM file |
| `ont-umi-lookup` | Match reads to UMI clusters from `ont-umi-cluster` output |
