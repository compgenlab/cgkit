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
- `MSA()` — incremental consensus multiple sequence alignment returning an `MSAAlignment` with optional homopolymer compression and reference sequence handling
- `MSAAlignment` — result type with `Consensus()`, `RehydratedConsensus()`, `WriteClustal()`, `WriteFasta()`, `GappedSequences()` for library-level output

### htsio — SAM/BAM/CRAM I/O

Native reading and writing of SAM, BAM, and tabix-indexed files. Samtools is only required for CRAM.

**Reading:**
- `SamReader` — interface with `Next()`, `Header()`, `Query()`, `Close()`
- `NewSamReader()` — auto-detects format: `.bam` → native BAM reader, `.sam`/`.sam.gz` → native text reader, `.cram` → samtools
- `Query(ref, start, end)` — returns `iter.Seq2[*SamRecord, error]` for indexed region queries (BAM via BAI, CRAM via samtools)
- Flag, MAPQ, and tag filtering via `SamReaderOpts`

**Writing:**
- `SamWriter` — interface with `Write()`, `Close()`
- `NewSamWriter()` — native BAM output (unsorted or coordinate/name sorted with merge sort), samtools for CRAM
- Sorted BAM writer buffers ~768MB, flushes to temp files, merge-sorts on Close

**Tabix:**
- `TabixReader` — query tabix-indexed BGZF files (BED, VCF, GFF) with TBI or CSI index auto-detection
- `TabixWriter` — sorted BGZF output with optional `.tbi` index generation; presets for BED, VCF, GFF
- Both use `iter.Seq2` for query results with 0-based half-open coordinates

**Index support:**
- BAI, TBI, CSI index parsers with shared `Query()` interface
- `ParseRegion()` — converts samtools-style region strings (`chr1:1000-2000`) to 0-based half-open

**Core types:**
- `SamRecord` — full SAM record with flag accessors (`IsUnmapped()`, `IsReverse()`, etc.) and typed tag access
- `SamHeader` — header manipulation including `@PG` line generation
- `TagFilter` — flexible tag-based filtering with comparison operators

### htsio/bgzf — BGZF compression

Low-level BGZF (Blocked GNU Zip Format) support used by BAM and tabix.

- `Reader` / `Writer` — streaming BGZF read/write with virtual offset tracking
- `IndexedReader` — random access with LRU block cache (default 64 blocks); supports virtual offset seeking and `.gzi` index for uncompressed offset seeking
- `NewBGZipFile()` — convenience constructor for file-backed BGZF output

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
| `seq-msa` | Multiple sequence alignment via incremental consensus (CLUSTAL by default; `--fasta` or `--consensus` for alternates; `--hp-compress` collapses homopolymers and rehydrates the consensus; `--ref <name>` marks a reference sequence that is aligned last, displayed first, and used for HP tiebreaks) |

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
