# cgltk

A toolkit for computational genomics, with a focus on Oxford Nanopore (ONT) long-read sequencing data.

---

## ont-umi-merge

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
