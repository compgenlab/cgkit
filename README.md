# cgltk

A toolkit for computational genomics, with a focus on Oxford Nanopore (ONT) long-read sequencing data.

---

## ont-umi-cluster

Collapses PCR and sequencing duplicates in a coordinate-sorted BAM file using UMI (Unique Molecular Identifier) barcodes. Reads with identical or near-identical UMIs originating from the same genomic locus are identified as duplicates; a consensus UMI sequence is computed for each group and written back to the BAM tags.

### Input

A coordinate-sorted BAM file with UMI sequences stored in a SAM tag (default `RX`). Two UMI formats are supported:

- Dash-separated: `XXXX-XXXX-XXXX-XXXX`
- TT-separated: `XXXXTTXXXXTTXXXXTTXXXX`

Both formats are treated as equivalent during comparison (see normalization below).

### Read grouping

In the default overlap mode, reads are grouped by genomic locus before clustering. Reads on the same strand whose alignments overlap or are within a configurable gap (`--overlap`, default 50 bp) are placed in the same group. UMI clustering is performed independently within each group. A whole-genome mode (`--whole-genome`) is also available, which treats all reads as a single group regardless of coordinate.

### UMI normalization

Before any comparison, UMI strings are normalized to dash-separated form by replacing each `TT` separator with `-`. This ensures that the two input formats are handled identically, and that a missed separator (e.g., `XXXX-XXXXXXXX-XXXX` instead of `XXXX-XXXX-XXXX-XXXX`) incurs an edit distance of 1 rather than being treated as an incomparable structure.

### Clustering algorithm

UMI clustering uses a graph-based single-linkage approach:

1. **Pairwise edit distances.** For all unique UMI strings within a group, the Levenshtein edit distance is computed for every pair. The edit distance operates on the full normalized (dash-separated) UMI string, so both base substitutions and single-base insertions/deletions are penalized uniformly. Structural variants such as a missed separator or an extra sequencing error in one group are therefore handled naturally without special-casing.

2. **Graph construction.** A graph is built with one node per unique UMI string. An undirected edge is added between two UMIs whose edit distance is at or below the threshold (`--umi-edit-distance`, default 3).

3. **Connected components.** Clusters are defined as the connected components of this graph, identified using a union-find (disjoint set union) data structure with path compression. Critically, two UMIs are placed in the same cluster if they are connected by *any path* of edges, not only if they are directly within the threshold of each other. This single-linkage criterion means that a chain A–B–C will be clustered together even if A and C are farther apart than the threshold, as long as each adjacent pair is within it.

   No assumption is made about which UMI is the "true" sequence based on read count. Count information does not influence cluster membership.

### Consensus calling

After clustering, a consensus UMI sequence is determined for each cluster by per-position majority vote:

- All member UMI strings are normalized to dash-separated form.
- Members whose normalized length matches the most common length in the cluster participate in the vote; length outliers are assigned to the cluster but do not influence the consensus sequence.
- At each position, the base with the highest total read count across all participating members is chosen. Votes are weighted by the number of reads carrying each UMI string, so higher-coverage UMIs have proportionally more influence on the consensus.

The resulting consensus sequence is always in dash-separated form.

### Output

- **BAM file** (`--output`, required): the input BAM with UMI tags updated. For each read whose UMI differs from its cluster consensus, the original UMI is moved to `--orig-umi-tag` (default `OX`) and the consensus is written to `--umi-tag` (default `RX`). Reads already matching the consensus are unchanged.

- **UMI counts TSV** (`--umi-counts`, optional): one row per unique UMI observed, with columns:
  - `chrom`, `start`, `end`, `strand` — genomic region of the overlap group (omitted in whole-genome mode)
  - `umi` — the original UMI string
  - `consensus` — the cluster consensus sequence
  - `count` — number of reads carrying this UMI
  - `consensus_total` — total reads in this cluster
  - `edit_dist` — Levenshtein distance from this UMI to the cluster consensus
  - `max_intra_cluster_dist` — maximum pairwise edit distance between any two members of this cluster; useful for evaluating cluster cohesion and tuning the edit distance threshold

- **BED file** (`--bed`, optional, overlap mode only): one row per overlap group processed, in BED6 format.

---

## seq-consensus-msa

Performs multiple sequence alignment of highly similar sequences using an incremental consensus algorithm. Designed for building consensus from ONT UMI read groups, but applicable to any set of closely related sequences.

### Input

A FASTA or FASTQ file containing the sequences to align. File format is detected by extension (`.fq`, `.fastq`, `.fq.gz`, `.fastq.gz` → FASTQ; otherwise FASTA). All sequences are loaded into memory. If the input is FASTQ, it is assumed that all reads in the file belong to the same alignment group.

### Algorithm

The alignment proceeds in three phases:

#### Phase 1: All-pairs pairwise alignment

Every pair of input sequences is aligned using global (Needleman-Wunsch) alignment. This produces N×(N-1)/2 alignments and their scores. The all-pairs step can be parallelized across multiple threads (`--threads`).

#### Phase 2: Seed pair selection

The pair of sequences with the **highest alignment score** is selected as the seed. Ties are broken by longest aligned length; if still tied, the first pair found is used. This pair forms the initial 2-sequence MSA. By selecting the best-scoring (and longest) pair, the initial alignment is anchored by the two most similar full-length reads.

#### Phase 3: Incremental incorporation

The remaining sequences are added one at a time:

1. Compute the **consensus** of the current MSA (majority vote at each column, ignoring gaps; ties broken alphabetically).
2. Align all unincorporated sequences to the consensus using **semi-global alignment** (the read is fully aligned end-to-end, but the consensus can have free end gaps on both sides — this accommodates truncated reads that don't span the full consensus).
3. Select the **best-scoring** unincorporated sequence and add it to the MSA.
4. Repeat until all sequences are incorporated.

Each addition refines the consensus, so each subsequent sequence is aligned to an increasingly accurate target. Because the highest-scoring sequences are added first, the consensus stabilizes quickly.

### Alignment modes

The command uses the Oxford Nanopore alignment preset by default (`--ont`, enabled by default), which includes:

- Match score: 1, mismatch penalty: 1
- Affine gap penalties tuned for ONT error profiles (insertion open: 2, deletion open: 3)
- Homopolymer indel discounts (reduced gap penalties in homopolymer runs)

An Illumina preset is available with `--ont=false`, which uses stricter mismatch and gap penalties appropriate for short-read error profiles.

Two alignment modes are used internally:

- **Global alignment** (Needleman-Wunsch) for the initial all-pairs distance computation — both sequences are fully aligned end-to-end.
- **Semi-global alignment** for the incorporation step — the query (read) is fully aligned, but the target (consensus) can have free end gaps. This handles truncated reads that are shorter than the consensus without incurring edge gap penalties.

### Consensus calling

The consensus sequence is built by majority vote at each column of the MSA:

1. For each column, count the occurrences of each non-gap base.
2. The base with the highest count is chosen as the consensus base.
3. Ties are broken alphabetically (A before C before G before T).
4. Columns where every entry is a gap are skipped.

The consensus is used in two places: internally during the incremental incorporation step (each new sequence is aligned to the current consensus), and as the final output when `--consensus` is specified.

The current implementation does not use FASTQ quality scores to weight the vote, nor does it produce quality scores for the consensus bases. Homopolymer-compressed consensus (tracking per-position run lengths and expanding back) is a planned future extension.

### Output

By default, the output is a **gapped multi-sequence FASTA** where each sequence includes `-` characters indicating gaps in the alignment. Sequence names are preserved from the input. Sequences appear in the order they were incorporated (seed pair first, then by descending alignment score).

With `--consensus`, the output is a single FASTA record containing the majority-vote consensus sequence.

### Usage

```
cgltk seq-consensus-msa <input.fasta|fastq> [flags]
```

### Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--ont` | `true` | Use Oxford Nanopore alignment defaults (set `--ont=false` for Illumina) |
| `-t`, `--threads` | `1` | Max parallel workers for the all-pairs alignment phase |
| `-o`, `--output` | stdout | Output file |
| `--consensus` | `false` | Output a single consensus sequence instead of the full MSA |
| `-v`, `--verbose` | `false` | Enable verbose alignment debug output |
