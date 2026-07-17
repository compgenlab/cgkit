# cgkit

CLI commands for computational genomics: sequence and alignment analysis, VCF
manipulation, BED/tabix wrangling, and NGS data operations, with particular focus
on Oxford Nanopore (long-read) sequencing workflows.

**Module:** `github.com/compgenlab/cgkit`

The underlying library (sequence I/O, alignment, SAM/BAM/CRAM handling) lives in
[`cghts`](https://github.com/compgenlab/cghts) (`github.com/compgenlab/cghts`).

## Building

```bash
make build     # Build all targets (darwin_arm64, darwin_amd64, linux_arm64, linux_amd64, windows_amd64)
make test      # Run all tests
```

Local development resolves the `cghts` dependency via a `go.work` workspace that
joins a sibling `cghts` checkout; release builds use the pinned module version in
`go.mod`.

## CLI commands

Usage: `cgkit [--profile=cpu.prof] <command>`

Run `cgkit <command> --help` for per-command flags. A `Since:` line in each
command's help shows the cgkit version it was added in.

### Oxford Nanopore

| Command | Description |
|---------|-------------|
| `ont-polya` | Find poly(A)/cleavage sites from a strand-specific aligned BAM |
| `ont-tags` | Find and trim common ONT tags from the start of reads in a FASTQ file |
| `ont-umi-cluster` | Collapse similar UMIs in a coordinate-sorted BAM file |
| `ont-umi-dedup` | Deduplicate UMI-clustered reads, keeping one representative per MI group |

### BED

| Command | Description |
|---------|-------------|
| `bed-clean` | Clean BED score entries to be integers (expands records to BED6+) |
| `bed-resize` | Resize BED regions (extend or shrink) |
| `bed-set` | Set algebra (intersect/union/subtract/exclusive) on two BED files |
| `bed-stats` | Summary statistics for a BED file |
| `bed-tobed3` | Convert a BED3+ file to a strict BED3 file |
| `bed-tobed6` | Convert a BED6+ file to a strict BED6 file |
| `bed-tofasta` | Extract FASTA sequences based on BED coordinates |

### FASTA/Q

| Command | Description |
|---------|-------------|
| `fasta-gc` | Return the GC content of sequences in a FASTA file |
| `fasta-wrap` | Reformat the sequences in a FASTA file to a specified line width |
| `fastq-gc` | Return the GC content of sequences in a FASTQ file |
| `fastq-tag` | Add a tag to the comment field of FASTQ records |

### SAM/BAM/CRAM

| Command | Description |
|---------|-------------|
| `sam-export` | Export columns and tags from a SAM/BAM/CRAM file as tab-delimited text |
| `sam-filter` | Filter SAM/BAM/CRAM reads and write to a new file |
| `sam-stats` | Summary statistics for a SAM/BAM/CRAM file |
| `sam-tofasta` | Convert SAM/BAM/CRAM reads to FASTA |
| `sam-tofastq` | Convert SAM/BAM/CRAM reads to FASTQ |

### Sequence

| Command | Description |
|---------|-------------|
| `seq-msa` | Multiple sequence alignment via incremental consensus |
| `seq-pairwise` | Align the two given sequences |
| `seq-revcomp` | Calculate the reverse-complement of the seq |

### Tabix

| Command | Description |
|---------|-------------|
| `tab-sort` | Sort a tab-delimited file and write as BGZF with optional tabix index |
| `tabix-index` | Build a tabix (.tbi) index for an existing BGZF-compressed file |

### VCF

| Command | Description |
|---------|-------------|
| `vcf-annotate` | Annotate a VCF file by adding INFO/FORMAT fields |
| `vcf-check` | Validate a VCF file |
| `vcf-chrfix` | Change the reference (chrom) format (Ensembl/UCSC) |
| `vcf-clearfilter` | Remove a filter from a VCF file |
| `vcf-concat` | Concatenate VCF files with the same samples but different variants |
| `vcf-export` | Export information from a VCF file as a tab-delimited file |
| `vcf-filter` | Filter a VCF file by stamping FILTER codes |
| `vcf-gtcount` | Summarize the genotype (GT) distribution across samples at given sites |
| `vcf-header-info` | Extract annotation/named fields from a VCF header |
| `vcf-merge` | Combine VCF files with the same variants but different annotations |
| `vcf-remove-flags` | Replace all INFO flags with a comma-separated list |
| `vcf-rename` | Change the names of samples |
| `vcf-reorder` | Reorder (or subset) the samples in a VCF file |
| `vcf-sample-export` | Write sample FORMAT values to a tab-delimited file, one sample per line |
| `vcf-samples` | Output the sample names in a VCF file |
| `vcf-split` | Split a VCF file into smaller files with N variants each |
| `vcf-stats` | Summary statistics about a VCF file |
| `vcf-strip` | Remove annotation and sample information, keeping VCF format |
| `vcf-svtofasta` | Extract SV breakend flanking sequences to FASTA |
| `vcf-tobed` | Export allele positions from a VCF file to BED format |
| `vcf-tobedpe` | Convert a structural-variant VCF to BEDPE format |
| `vcf-tocount` | Convert a VCF to a count file using the AD (or RO/AO) format field |
| `vcf-tstv` | Calculate a Ts/Tv ratio for SNVs |
