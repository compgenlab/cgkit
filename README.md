# cgkit

CLI commands for computational genomics: sequence and alignment analysis, VCF
manipulation, BED/tabix wrangling, and NGS data operations.

Oxford Nanopore UMI and poly(A) tooling lives in
[`nupa`](https://github.com/compgenlab/nupa), a focused toolkit built on the same
library; the `ont-*` commands that used to be here moved there.

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

## Remote inputs

Inputs may be local paths, `http(s)://` URLs, or `s3://` objects:

```bash
cgkit sam-stats s3://bucket/sample.bam
cgkit vcf-varquery https://host/cohort.vcf.gz --variant chr1:1000:A:T
cgkit vcf-varquery s3://bucket/cohort/ --sample NA12878
```

Indexes (`.bai`, `.crai`, `.tbi`, `.csi`, `.fai`) are resolved over the same
transport as the data, so an indexed region query still seeks and transfers only
the byte ranges it needs. `--cram-ref` takes a locator too, independently of
where the reads live.

S3 credentials come from the standard AWS chain — environment, shared config
(honouring `AWS_PROFILE`), then an instance or container role; `aws sts
get-caller-identity` is the check. `AWS_ENDPOINT_URL` points at an S3-compatible
gateway such as MinIO or Ceph.

Two costs worth knowing. A *streaming* read of a remote object transfers the
whole thing, because a stream has no index to skip with — commands taking
`--region`, and any query against a Parquet store, seek instead. And an
unindexed remote VCF cannot seek at all, so every query re-reads the object;
`vcf-varsummary` reports whether an index was found.

**Outputs are always local.** A remote `-o` or `--out` is refused by name rather
than failing inside a writer.

## CLI commands

Usage: `cgkit [--profile=cpu.prof] <command>`

`--profile=FILE` writes a CPU profile and may appear anywhere on the command
line, before or after the subcommand. It is stripped from the recorded
invocation, so the provenance a run stamps into its output describes a
reproducible command rather than a profiling one.

Run `cgkit <command> --help` for per-command flags. A `Since:` line in each
command's help shows the cgkit version it was added in.

A command invoked with no arguments is a usage error on stderr with a non-zero
exit, not a help page on stdout — so `cgkit vcf-tobed > out.bed && next-step
out.bed` cannot proceed with help text as its input.

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
| `vcf-toparquet` | Convert a VCF to a sparse Parquet genotype store ([format docs](docs/vcf-toparquet.md)) |
| `vcf-tstv` | Calculate a Ts/Tv ratio for SNVs |
| `vcf-varquery` | Query genotypes by site, by sample, or both ([format docs](docs/vcf-toparquet.md)) |
| `vcf-varsummary` | Describe a store or VCF: samples, contigs, provenance, per-chromosome census |
