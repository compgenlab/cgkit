# cgkit

CLI commands for computational genomics: sequence analysis, NGS data wrangling,
and bioinformatics operations, with particular focus on Oxford Nanopore
(long-read) sequencing workflows.

**Module:** `github.com/compgenlab/cgkit`

The underlying library (sequence I/O, alignment, SAM/BAM/CRAM handling) lives in
[`hts`](https://github.com/compgenlab/hts) (`github.com/compgenlab/hts`).

## Building

```bash
make build     # Build all targets (darwin_arm64, linux_arm64, linux_amd64)
make test      # Run all tests
```

Local development resolves the `hts` dependency via a `go.work` workspace that
joins a sibling `hts` checkout; release builds use the pinned module version in
`go.mod`.

## CLI commands

Usage: `cgkit [--profile=cpu.prof] <command>`

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
| `sam-tofasta` | Convert reads to FASTA (optionally writing SAM tags into the comment) |
| `sam-tofastq` | Convert reads to FASTQ (optionally writing SAM tags into the comment) |

### Oxford Nanopore

| Command | Description |
|---------|-------------|
| `ont-polya` | Find per-read poly(A)/cleavage sites in a strand-specific aligned BAM |
| `ont-tags` | Find and trim common ONT adapter/primer tags from FASTQ reads |
| `ont-umi-cluster` | Collapse similar UMIs in a coordinate-sorted BAM file |
| `ont-umi-dedup` | Keep one representative read per UMI group (`MI` tag) from a coordinate-sorted BAM |
| `ont-umi-lookup` | Match reads to UMI clusters from `ont-umi-cluster` output |
