# ont-umi-cluster

Collapses PCR and sequencing duplicates in a coordinate-sorted BAM file using UMI (Unique Molecular Identifier) barcodes. Reads with identical or near-identical UMIs originating from the same genomic locus are clustered together; a representative UMI sequence is selected for each cluster and written back to the BAM tags.

## Usage

```
cgltk ont-umi-cluster [flags] <input.bam>
```

Required: `--output` (`-o`) specifying the output BAM file path.

## Input

A coordinate-sorted BAM file with UMI sequences stored in a SAM tag (default `RX`). Three UMI separator formats are supported:

- Slash-separated: `XXXX/XXXX/XXXX/XXXX`
- Dash-separated: `XXXX-XXXX-XXXX-XXXX`
- TT-separated: `XXXXTTXXXXTTXXXXTTXXXX`

All formats are normalized to `/`-separated form before comparison.

## Pipeline overview

The pipeline has three phases:

1. **Read grouping** (detection): stream reads from the BAM, group overlapping reads into components using a union-find with bin-indexed overlap detection.
2. **UMI clustering**: within each component, compute all-pairs edit distances between unique UMIs and cluster them using the selected method.
3. **Output**: write the BAM with updated UMI tags and optional per-cluster counts.

Phases 1 and 2 are decoupled via a worker pool: the main goroutine handles detection (single-threaded, I/O-bound) and completed components are sent to N worker goroutines for parallel UMI clustering (CPU-bound).

## Read grouping

Reads are grouped by genomic position using bin-indexed overlap detection. Two reads are candidates for the same group if they satisfy the overlap criteria on the same strand (or any strand with `--no-strand`).

### Default mode (match-both-ends)

Two reads match if their 3' ends are within `--overlap` bp of each other. The 5' proximity is enforced implicitly by the coordinate-sorted ejection: reads whose start position falls more than `--overlap` behind the current position are ejected from the active set.

### `--match-one-end` mode

Two reads match if EITHER their 5' starts OR their 3' ends are within `--overlap` bp. Reads are ejected when the current position passes their 3' end by more than `--overlap`. This mode is used for samples with 5' degradation where reads from the same molecule may not share a common 5' start.

### Bin-indexed detection

Instead of a linear scan of the active buffer for each incoming read (O(buffer_size) per read), the active reads are indexed in two bin maps keyed by `position / overlap`:

- `endIndex`: bins reads by 3' end position
- `startIndex`: bins reads by 5' start position

Each incoming read queries at most 3 bins (target +/- 1) per index. Combined with a root-skip optimization (once a mega-component has formed, bin entries whose union-find root matches the current read's root are skipped after a single `find()` call), detection is O(matches) per read instead of O(buffer_size).

### `--region` mode

Process only a single samtools region (e.g., `chr19` or `chr19:1000-2000`). Skipped-ref and unmapped read pass-through is disabled in this mode, allowing the user to orchestrate per-region jobs externally (e.g., one SLURM task per chromosome). Without `--region`, the entire BAM is read in a single pass with automatic chromosome-transition detection.

## UMI clustering

### Edit distance

UMI similarity is measured by Levenshtein edit distance on the normalized (slash-separated) UMI string. The implementation uses a bounded DP (Ukkonen's cutoff): after filling each row, if the row minimum exceeds the threshold, the computation bails out immediately. For a threshold of 3 on 19-character UMIs, most dissimilar pairs exit after 3-4 DP rows instead of the full 19x19 matrix.

### HP-aware edit distance (`--hp-dist`)

When `--hp-dist` is set, homopolymer indels cost 0 instead of 1. Specifically:

- Deleting `a[i]` costs 0 if `a[i] == a[i-1]` (shrinking an HP run)
- Inserting `b[j]` costs 0 if `b[j] == b[j-1]` (shrinking an HP run)
- Substitutions always cost 1

This correctly handles ONT's most common error mode (homopolymer-length variation) without the distortion that full homopolymer compression causes. For example:

| UMI_A | UMI_B | Standard distance | HP-aware distance |
|-------|-------|-------------------|-------------------|
| AAGA  | AAAA  | 1 (G to A)       | 1 (substitution preserved) |
| AAACG | AACG  | 1 (HP indel)     | 0 (HP indel discounted) |
| AAGA  | AAGGA | 1 (HP indel)     | 0 (HP indel discounted) |

Full HP compression would distort the first case: `AAGA` compresses to `AGA`, `AAAA` to `A`, inflating the distance from 1 to 2. The HP-aware distance avoids this by modifying the DP cost model rather than the input strings.

### Clustering methods (`--umi-cluster-method`)

All methods operate on the same set of edges (pairs of UMIs within the edit distance threshold). They differ in how edges are used to form clusters.

#### `connected` (single-linkage)

Union-find on all edges. Two UMIs are in the same cluster if they are connected by any path of edges, regardless of path length. This is the classical approach but suffers from **chaining**: A-B-C can form one cluster even if d(A,C) >> threshold.

#### `adjacency` (default)

Greedy assignment. UMIs are processed in count-descending order. Each unassigned UMI becomes a cluster center; all its unassigned direct neighbors (within threshold) join that cluster. No chaining occurs because each UMI is assigned exactly once and merges are one-hop only.

This is the recommended method for ONT data. Two high-count UMIs that are close (e.g., 2 edits apart) always merge — whichever has the higher count grabs the other. This handles the ONT error model where two well-amplified copies of the same molecule can have independent sequencing errors.

#### `directional`

Edges are filtered by a PCR error count model before union-find:

```
count(low) <= 2 * count(high) * (1/4)^distance
```

Only low-count UMIs that could plausibly be PCR/sequencing errors of a high-count UMI create edges. Two equally-expressed UMIs do not merge. This model is from UMI-tools (Smith et al., Genome Research 2017) and is calibrated for Illumina error rates (~0.1%). It is generally **too conservative for ONT** where per-base error rates are 5-15%.

#### `tiered` (distance-attenuated BFS clustering)

A novel method that combines the center-selection from adjacency with multi-hop expansion at decreasing stringency. Starting from the highest-count unassigned UMI (the center), a BFS expands outward with the allowed edit distance decreasing by 1 at each hop:

```
With --umi-edit-distance 3:
  Hop 0: center (highest-count unassigned UMI)
  Hop 1: merge neighbors at d <= 3
  Hop 2: merge neighbors at d <= 2
  Hop 3: merge neighbors at d <= 1
  Stop.
```

This is more permissive than adjacency (which is one hop only) but prevents chaining at high distances: a chain center-A(d=3)-B(d=3) is blocked because hop 2 requires d <= 2. Very similar UMIs (d=1) can chain up to 3 hops deep, capturing errors-of-errors, while distant UMIs (d=3) are limited to direct neighbors of the center.

### Adaptive threshold (`--adaptive-threshold`)

A statistical post-filter that discards edges at distances where random collisions are likely to dominate.

#### Rationale

When clustering UMIs, the edit distance threshold determines which pairs of UMIs are considered "similar enough" to potentially originate from the same molecule. However, the number of random collisions (pairs of independent UMIs that happen to be within edit distance d by chance) grows with both the component size and the distance:

- In a small component (100 UMIs), a pair at distance 3 is almost certainly real — the probability of a random collision is negligible.
- In a large component (100,000 UMIs), there are ~5 billion possible pairs. Even at the low per-pair collision probability of ~3.5 x 10^-6 for distance 3, the expected number of random close pairs is ~17,600. These false edges, when fed into any clustering method, cause incorrect merges.

Rather than choosing a fixed threshold that works for all component sizes, the adaptive threshold lets the data decide: at each edit distance, it measures how many edges were actually found and compares against the expected number of random collisions. If the false positive rate exceeds a tolerance (default 5%), that distance is excluded.

#### Method

After all-pairs edge-finding, the adaptive filter computes a per-distance false positive rate. The expected number of random pairs at exactly distance d between N independent UMIs of length L over a 4-letter alphabet is:

```
P(exactly d) = C(L, d) * 3^d / 4^L
E_false(d) = N * (N-1) / 2 * P(exactly d)
```

Where:
- `C(L, d)` is "L choose d" — the number of ways to pick d positions to mutate
- `3^d` accounts for the 3 possible wrong bases at each mutated position
- `4^L` is the total UMI sequence space

This is a substitution-only approximation. Indels slightly expand the collision neighborhood, so the true collision probability is marginally higher — making this a conservative (slightly permissive) estimate.

The false positive rate at distance d is then:

```
FPR(d) = E_false(d) / actual_edges(d)
```

where `actual_edges(d)` is the number of edges found at exactly distance d by the all-pairs computation. If `FPR(d)` exceeds `--adaptive-alpha` (default 0.05), all edges at distance d are discarded before clustering.

#### Concrete example (L=16 bases)

The per-pair collision probabilities are:

| Distance d | P(exactly d) | Interpretation |
|---|---|---|
| 1 | 1.1 x 10^-8 | ~1 in 90 million |
| 2 | 2.5 x 10^-7 | ~1 in 4 million |
| 3 | 3.5 x 10^-6 | ~1 in 286,000 |

Expected false pairs at each component size:

| N (unique UMIs) | d=1 expected false | d=2 expected false | d=3 expected false |
|---|---|---|---|
| 100 | ~0 | ~0 | 0.02 |
| 1,000 | 0.006 | 0.13 | 1.8 |
| 10,000 | 0.6 | 13 | 176 |
| 100,000 | 56 | 1,257 | 17,600 |

At 5% FPR, the minimum number of real edges needed at each distance to survive filtering is `E_false / 0.05`:

| N | d=1 min edges | d=2 min edges | d=3 min edges |
|---|---|---|---|
| 100 | 1 | 1 | 1 |
| 1,000 | 1 | 3 | 36 |
| 10,000 | 12 | 252 | 3,520 |
| 100,000 | 1,118 | 25,140 | 352,000 |

In practice: components with N <= 1,000 keep all distances. Components with N >= 10,000 progressively lose higher distances. A component with 86,000 unique UMIs (observed in real ONT data) had 17,600 expected random pairs at d=3 — these false edges caused single-linkage to create an 81,655-member mega-cluster. The adaptive threshold would have excluded d=3 (and possibly d=2) edges, preventing the false chaining entirely.

#### Properties

- **Independent of clustering method**: the adaptive threshold is a pre-filter on edges, applied before any of the four clustering methods (connected, adjacency, directional, tiered).
- **Transparent**: when edges are excluded, a diagnostic message is printed to stderr showing the distance, FPR, edge count, and expected false count.
- **Conservative default**: the 5% FPR default means up to 1 in 20 edges at a given distance may be a random collision. This is permissive enough to avoid discarding real edges in small components, while aggressive enough to catch the random-collision problem in large ones.
- **Tunable**: `--adaptive-alpha` controls the tolerance. Lower values (e.g., 0.01) are more conservative (discard more aggressively); higher values (e.g., 0.10) are more permissive.

### Representative selection

The representative UMI for each cluster is the highest-count member (ties broken by UMI string length, then lexicographic order). The representative is always written in normalized (`/`-separated) form.

## Output

### BAM file (`--output`, required)

The input BAM with UMI tags updated. For each read whose UMI differs from its cluster's representative:
- The original UMI is moved to `--tag-orig` (default `OX`)
- The representative is written to `--tag-umi` (default `RX`)

Reads already matching the representative are unchanged. The output is coordinate-sorted via samtools sort.

### Molecule ID (`--tag-mi`)

When `--tag-mi` is set, each read receives an `MI` tag with a unique molecule group identifier (e.g., `mi_000000001`). All reads in the same UMI cluster share the same MI value.

### UMI counts (`--summary-counts`, optional)

A tab-delimited BED6+ file with one row per UMI cluster per read-overlap-group. Coordinates are per-cluster (the bounding box of reads carrying that cluster's UMIs), not per-component. This file is not sorted by default. Columns:

| # | Column | Description |
|---|--------|-------------|
| 1 | chrom | Reference name |
| 2 | start | 0-based start of the cluster's reads |
| 3 | end | 0-based end (exclusive) of the cluster's reads |
| 4 | name | Molecule ID (e.g., `mi_000000001`) |
| 5 | score | Total read count in the cluster |
| 6 | strand | `+`, `-`, or `.` |
| 7 | representative | Chosen representative UMI (slash-separated) |
| 8 | numUMIs | Number of distinct original UMIs in the cluster |
| 9 | maxEditDist | Max pairwise edit distance within the cluster (-1 if skipped for large clusters) |
| 10 | umis | Comma-separated list of original UMIs |

## Performance

### Parallelism

The `--threads` (`-t`) flag controls the number of worker goroutines for UMI clustering. The main detection loop runs on a single goroutine; completed components are dispatched to workers via a buffered channel. Within each worker, the all-pairs edge-finding and the max-intra-cluster-distance computation are both parallelized across the same thread budget using round-robin row distribution.

samtools sort (the output writer) uses 2 fixed threads and does not consume significant CPU.

### Bounded edit distance

The Levenshtein DP uses Ukkonen's cutoff: after each row, if the row minimum exceeds the bound, the function returns immediately. For a threshold of 3 on 19-character UMIs, the vast majority of pairs bail after 3-4 rows. This provides a 5-10x speedup over unbounded DP.

### Max-intra-cluster-distance cap

The all-pairs max-intra-cluster-distance computation is capped at 10,000 members. Clusters exceeding this size report `-1` instead of computing the O(n^2) metric, which would take hours for clusters with 80,000+ members.

## Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--adaptive-alpha` | 0.05 | Maximum false positive rate per edit distance level |
| `--adaptive-threshold` | false | Discard edges at distances exceeding the false positive rate threshold |
| `--hp-dist` | false | HP-aware edit distance (HP indels cost 0) |
| `--ignore-refs` | | References to pass through without clustering (comma-separated) |
| `--match-one-end` | false | Match reads if EITHER 5' or 3' ends are within gap |
| `--no-strand` | false | Ignore strand when grouping reads |
| `-o`, `--output` | (required) | Output BAM file path |
| `--overlap` | 50 | Maximum gap (bp) between reads to group them |
| `--region` | | Process only this region |
| `--tag-mi` | false | Add MI tag with molecule group ID |
| `--tag-orig` | `OX` | SAM tag to store original UMI before correction |
| `--tag-umi` | `RX` | SAM tag containing UMI sequence |
| `-t`, `--threads` | 1 | Worker threads for UMI clustering |
| `--summary-counts` | | Write per-cluster UMI summary to this file |
| `--umi-cluster-method` | `adjacency` | Clustering method: `connected`,`adjacency`, `directional`, `tiered` |
| `--umi-edit-distance` | 3 | Maximum Levenshtein edit distance to cluster two UMIs |
`--whole-genome` | false | Process all UMIs as a single group (ignore coordinates) |

## References

- Smith T, Heger A, Sudbery I. "UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy." *Genome Research* 27(3):491-499, 2017. (Adjacency and directional clustering methods.)
