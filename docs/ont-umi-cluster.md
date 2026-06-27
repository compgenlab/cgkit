# ont-umi-cluster

Collapses PCR and sequencing duplicates in a coordinate-sorted BAM file using UMI (Unique Molecular Identifier) barcodes. Reads with identical or near-identical UMIs originating from the same genomic locus are clustered together; a representative UMI sequence is selected for each cluster and written back to the BAM tags.

## Usage

```
cgkit ont-umi-cluster [flags] <input.bam>
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

### Splice junction matching (`--junction-match`)

When `--junction-match` is set, reads must have compatible splice junctions (from CIGAR `N` operations) in addition to positional overlap to be grouped together. This is intended for RNA-seq data where reads from different transcript isoforms may overlap positionally but represent distinct molecules.

#### Junction extraction and pre-merging

Splice junctions are extracted from each read's CIGAR string. Each `N` operation produces a junction with a donor position (where the intron starts) and an acceptor position (where the next exon begins). Adjacent junctions separated by ≤ `--junction-window` bp (default 20) are pre-merged into a single spanning junction. This handles cases where a very small exon is present in one read's alignment but missed in another — a common occurrence with Oxford Nanopore alignments. For example:

```
Read 1: 50M──N(1000bp)──8M──N(500bp)──50M   → 2 junctions
Read 2: 50M──N(1508bp)──50M                  → 1 junction
```

After pre-merging (gap of 8bp ≤ 10bp window), read 1's two junctions collapse into one spanning junction that matches read 2's single junction.

#### Compatibility rules

- **No junctions + no junctions**: compatible (unspliced reads group together).
- **No junctions + has junctions**: incompatible (an unspliced read and a spliced read at the same locus are likely different molecules).
- **Both have junctions, default mode**: junction sets must match exactly — same count, with each junction's donor and acceptor positions within ±`--junction-window` bp.
- **Both have junctions, `--match-one-end` mode**: one read's junction set may be a contiguous sub-sequence of the other's, anchored at the matching end. When 3' ends match (suffix anchor), the shorter read's junctions must match a suffix of the longer read's junctions. When 5' starts match (prefix anchor), the shorter read's junctions must match a prefix. This handles 5' truncation where the shorter molecule is missing junctions from the truncated end.

#### Example (`--match-one-end` with 5' truncation)

```
Read 1 (full):      junctions A, B, C, D     start=100, end=5000
Read 2 (truncated): junctions    B, C, D     start=800, end=5000
```

The 3' ends match (both at 5000). Read 2's junctions [B, C, D] are a suffix of read 1's [A, B, C, D] — compatible. Junction A is missing because read 2 starts after the first exon-intron boundary.

## UMI clustering

### Edit distance

UMI similarity is measured by Levenshtein edit distance on the normalized (slash-separated) UMI string. The implementation uses a bounded DP (Ukkonen's cutoff): after filling each row, if the row minimum exceeds the threshold, the computation bails out immediately. For a threshold of 3 on 19-character UMIs, most dissimilar pairs exit after 3-4 DP rows instead of the full 19x19 matrix.

### HP-aware edit distance (`--hp-dist`)

When `--hp-dist` is set, one HP indel per UMI segment is free. The free is shared between insertions and deletions — only one HP indel per segment is discounted, regardless of direction. The free resets when either string crosses a `/` separator boundary. Substitutions and non-HP indels always cost 1.

This is implemented as an augmented Levenshtein DP with state `dp[i][j][f]` where `f ∈ {0, 1}` tracks whether the single free HP indel for the current segment has been consumed.

This correctly handles ONT's most common error mode (±1 homopolymer-length variation) while preventing long HP runs from collapsing arbitrarily. For example:

| UMI_A | UMI_B | Standard distance | HP-aware distance | Notes |
|-------|-------|-------------------|-------------------|-------|
| AAGA  | AAAA  | 1 | 1 | Substitution preserved |
| AAACG | AACG  | 1 | 0 | ±1 in A-run, free |
| AAGA  | AAGGA | 1 | 0 | ±1 in G-run, free |
| AACC  | AAC   | 1 | 0 | ±1 in A-run, free |
| AACC  | AC    | 2 | 1 | One free + one paid (shared budget per segment) |
| CCGA  | CGAA  | 2 | 1 | One free (del C) + one paid (ins A), shared |
| AAAA  | AA    | 2 | 1 | ±1 free, second HP indel costs 1 |
| AAAA  | A     | 3 | 2 | ±1 free, two more cost 1 each |
| GGGG  | G     | 3 | 2 | Long HP collapse no longer free |
| CGGC  | CCCC  | 2 | 2 | Pure substitutions, no HP discount |

Full HP compression distorts substitutions: `AAGA` compresses to `AGA`, `AAAA` to `A`, inflating the distance from 1 to 2. Making all HP indels cost 0 (as in an earlier version of `--hp-dist`) is too lenient — sequences like `AAAA/GGGG/CCGA/CCCC` and `CCCC/AAAA/CCGA/GGGG` (standard distance 12) appeared at HP-distance 3 because all HP indels were free, causing false clustering. The per-segment first-free model avoids both problems: substitutions are preserved, ±1 HP errors are tolerated, and long HP collapses are penalized.

### Clustering methods (`--umi-cluster-method`)

All methods operate on the same set of edges (pairs of UMIs within the edit distance threshold). They differ in how edges are used to form clusters.

#### `connected` (single-linkage)

Union-find on all edges. Two UMIs are in the same cluster if they are connected by any path of edges, regardless of path length. This is the classical approach but suffers from **chaining**: A-B-C can form one cluster even if d(A,C) >> threshold.

#### `adjacency` (default)

Greedy assignment. UMIs are processed in count-descending order. Each unassigned UMI becomes a cluster center; all its unassigned direct neighbors (within threshold) join that cluster. No chaining occurs because each UMI is assigned exactly once and merges are one-hop only.

This is the recommended method for ONT data. Two high-count UMIs that are close (e.g., 2 edits apart) always merge — whichever has the higher count grabs the other. This handles the ONT error model where two well-amplified copies of the same molecule can have independent sequencing errors.

Because adjacency is strictly one-hop, UMIs that are only reachable through intermediate nodes are not clustered, even if those intermediates are close. For example, with `--umi-edit-distance 3` and these edges:

```
A(100)—d=1—B(40)    B—d=1—C(15)    C—d=2—D(8)    A—d=3—D    D—d=1—E(3)
```

A becomes center and grabs its direct neighbors B (d=1) and D (d=3). C and E have no edge to A, so they are not claimed. C becomes its own center next (count 15), but its neighbors B and D are already assigned — so C is a singleton. E likewise becomes a singleton.

Result: {A, B, D}, {C}, {E}.

By contrast, tiered with the same edges produces {A, B, C, D, E}: B and D each pull in their neighbors on hop 2 (with `maxDist = 2`), reaching C and E that adjacency cannot.

#### `directional`

Edges are filtered by a PCR error count model before union-find:

```
count(low) <= 2 * count(high) * (1/4)^distance
```

Only low-count UMIs that could plausibly be PCR/sequencing errors of a high-count UMI create edges. Two equally-expressed UMIs do not merge. This is a two-step process: first, edges that fail the count-ratio test are removed; then, union-find is applied to the surviving edges with full transitive closure (the same as `connected`). Chaining can still occur — if A–B and B–C both pass the count filter, A and C end up in the same cluster regardless of d(A,C) — but the count filter makes it harder for edges to survive, reducing chaining in practice.

This model is from UMI-tools (Smith et al., Genome Research 2017) and is calibrated for Illumina error rates (~0.1%). It is generally **too conservative for ONT** where per-base error rates are 5-15%.

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

##### Worked example

Consider two independent groups of UMIs (e.g., from different molecules) with `--umi-edit-distance 3`:

```
Group 1                          Group 2
UMI  Count                       UMI  Count
 A    100                         L    90
 B     40                         M    35
 C     15                         N    12
 D      8                         O     6
 E      3                         P     2

Edges (group 1):                 Edges (group 2):
 A—B: d=1                        L—M: d=1
 B—C: d=1                        M—N: d=1
 C—D: d=2                        N—O: d=2
 A—D: d=3                        L—O: d=3
 D—E: d=1                        O—P: d=1
```

No edges exist between the two groups.

**Processing order** (by count descending): A(100), L(90), B(40), M(35), ...

**Step 1 — A becomes center, BFS runs to completion:**

| Queue entry | maxDist (= 3 − hop) | Outcome |
|---|---|---|
| A, hop 0 | 3 | Neighbors: B (d=1 ≤ 3 ✓), D (d=3 ≤ 3 ✓). Both join cluster A at hop 1. |
| B, hop 1 | 2 | Neighbors: C (d=1 ≤ 2 ✓). C joins at hop 2. |
| D, hop 1 | 2 | Neighbors: C (assigned, skip), E (d=1 ≤ 2 ✓). E joins at hop 2. |
| C, hop 2 | 1 | Neighbors: D (assigned, skip). No new members. |
| E, hop 2 | 1 | No unassigned neighbors. |

Cluster A = {A, B, C, D, E}. BFS is complete.

**Step 2 — L becomes center, BFS runs to completion:**

The same process repeats independently: L claims M, N, O, P via BFS with the same hop logic. Cluster L = {L, M, N, O, P}.

**Result:** Two separate clusters. No merging occurs between them — clusters are never combined after formation.

##### Comparison to other methods on a chain topology

Now consider a chain where every link is d=3:

```
A (100) —d=3— B (5) —d=3— C (2)
```

| Method | Result | Why |
|---|---|---|
| **connected** | {A, B, C} | Union-find merges transitively: A–B and B–C ⇒ A–B–C. |
| **adjacency** | {A, B}, {C} | One hop only: A grabs B, but C is not A's neighbor. |
| **tiered** | {A, B}, {C} | A→B at hop 1 uses full budget (d=3 ≤ 3). B→C at hop 2 requires d ≤ 2, but d=3 > 2. Blocked. |

If the chain uses d=1 edges instead:

```
A (100) —d=1— B (5) —d=1— C (2)
```

| Method | Result | Why |
|---|---|---|
| **connected** | {A, B, C} | Same as before. |
| **adjacency** | {A, B}, {C} | Still one hop — C is not a direct neighbor of A. |
| **tiered** | {A, B, C} | A→B at hop 1 (d=1 ≤ 3). B→C at hop 2 (d=1 ≤ 2). Allowed. |

This illustrates the key property of tiered: close UMIs (d=1) can chain through multiple hops, capturing errors-of-errors, while distant UMIs (d=3) cannot chain beyond the center.

### Adaptive threshold (`--adaptive-threshold`)

A statistical post-filter that discards edges at distances where random collisions are likely to dominate.

#### Rationale

When clustering UMIs, the edit distance threshold determines which pairs of UMIs are considered "similar enough" to potentially originate from the same molecule. However, the number of random collisions (pairs of independent UMIs that happen to be within edit distance d by chance) grows with both the component size and the distance:

- In a small component (100 UMIs), a pair at distance 3 is almost certainly real — the probability of a random collision is negligible.
- In a large component (100,000 UMIs), there are ~5 billion possible pairs. Even at the low per-pair collision probability of ~3.5 x 10^-6 for distance 3, the expected number of random close pairs is ~17,600. These false edges, when fed into any clustering method, cause incorrect merges.

Rather than choosing a fixed threshold that works for all component sizes, the adaptive threshold lets the data decide: at each edit distance, it measures how many cumulative edges were actually found and compares against the expected number of random collisions. If the cumulative false positive rate exceeds a tolerance (default 10%), that distance and all higher distances are excluded.

#### Method

After all-pairs edge-finding, the adaptive filter computes a cumulative false positive rate. The cumulative probability that two independent random UMIs of length L over a 3-letter alphabet (as used by ONT UMIs) are within edit distance d is:

```
P(≤d) = Σ_{k=0}^{d} C(L, k) * 2^k / 3^L
```

Where:
- `C(L, k)` is "L choose k" — the number of ways to pick k positions to mutate
- `2^k` accounts for the 2 possible wrong bases at each mutated position (3-letter alphabet)
- `3^L` is the total UMI sequence space

This is a substitution-only approximation. Indels slightly expand the collision neighborhood, so the true collision probability is marginally higher — making this a conservative (slightly permissive) estimate.

The expected number of false edges at cumulative distance ≤d among N UMIs is:

```
E_false(≤d) = N * (N-1) / 2 * P(≤d)
```

The cumulative false positive rate is then:

```
FPR(≤d) = E_false(≤d) / actual_edges(≤d)
```

where `actual_edges(≤d)` is the total number of edges found at distances 1 through d. If `FPR(≤d)` exceeds `--adaptive-alpha` (default 0.10), all edges at distance d and above are discarded before clustering.

Because the filter uses cumulative counts, it is monotonic: once the FPR exceeds the threshold at some distance d, all higher distances are automatically excluded. This prevents the anomaly where a per-distance filter might exclude d=2 but keep d=3 due to an unusually high number of real edges at d=3.

#### Concrete example (L=16 bases, 3-letter alphabet)

The cumulative collision probabilities are:

| Distance | P(≤d) | Interpretation |
|---|---|---|
| ≤1 | 7.7 x 10^-7 | ~1 in 1.3 million |
| ≤2 | 1.2 x 10^-5 | ~1 in 84,000 |
| ≤3 | 1.2 x 10^-4 | ~1 in 8,600 |

Each level is ~10-15× more likely than the previous.

Expected cumulative false edges at each component size:

| N (unique UMIs) | nPairs | E_false(≤1) | E_false(≤2) | E_false(≤3) |
|---|---|---|---|---|
| 100 | 4,950 | 0.004 | 0.06 | 0.6 |
| 1,000 | 499,500 | 0.4 | 6 | 58 |
| 5,000 | 12.5M | 10 | 149 | 1,450 |
| 10,000 | 50M | 38 | 595 | 5,800 |
| 20,000 | 200M | 153 | 2,380 | 23,200 |

At 10% FPR (default alpha), the minimum number of cumulative edges needed to survive filtering is `E_false / 0.10`:

| N | min edges for d≤1 | min edges for d≤2 | min edges for d≤3 |
|---|---|---|---|
| 100 | 1 | 1 | 6 |
| 1,000 | 4 | 59 | 580 |
| 5,000 | 96 | 1,487 | 14,500 |
| 10,000 | 383 | 5,950 | 58,000 |
| 20,000 | 1,533 | 23,800 | 232,000 |

The 3-letter alphabet (ONT UMIs use only A, C, G) has a much smaller sequence space (3^16 ≈ 43 million) compared to a 4-letter alphabet (4^16 ≈ 4.3 billion), making random collisions 20-67× more likely. With the default alpha of 0.10 and ~1-2 real neighbors per UMI, the adaptive threshold begins excluding d≤3 edges at around N ≈ 2,000-3,000 unique UMIs.

#### Worked example

Consider a component with 100 unique UMIs (N=100) and UMI length L=12, with `--umi-edit-distance 3`:

```
nPairs = 100 × 99 / 2 = 4,950
```

Suppose the all-pairs computation found these edges:

| Distance | Edges at this distance | Cumulative edges |
|---|---|---|
| 1 | 15 | 15 |
| 2 | 8 | 23 |
| 3 | 5 | 28 |

The cumulative FPR at each level:

| Level | E_false(≤d) | Cumulative edges | FPR |
|---|---|---|---|
| ≤1 | 0.22 | 15 | 0.22/15 = **1.5%** |
| ≤2 | 2.4 | 23 | 2.4/23 = **10.4%** |
| ≤3 | 15.3 | 28 | 15.3/28 = **54.7%** |

With alpha=10%:
- **d≤1**: FPR 1.5% < 10% — keep
- **d≤2**: FPR 10.4% > 10% — **exclude d=2 and d=3**

The effective threshold drops to 1. All d=2 and d=3 edges are discarded before clustering. Because the filter is cumulative, there is no possibility of excluding d=2 while keeping d=3.

#### Properties

- **Independent of clustering method**: the adaptive threshold is a pre-filter on edges, applied before any of the four clustering methods (connected, adjacency, directional, tiered).
- **Monotonic**: uses cumulative FPR, so once a distance is excluded, all higher distances are excluded too.
- **Transparent**: when edges are excluded, a diagnostic message is printed to stderr showing the distance, cumulative FPR, cumulative edge count, and expected false count.
- **Permissive default**: the 10% FPR default means up to 1 in 10 cumulative edges may be a random collision. This is permissive enough to avoid discarding real edges in small components, while aggressive enough to catch the random-collision problem in large ones.
- **Tunable**: `--adaptive-alpha` controls the tolerance. Lower values (e.g., 0.01) are more conservative (discard more aggressively); higher values (e.g., 0.20) are more permissive.

### Representative selection

The representative UMI for each cluster is the highest-count member (ties broken by UMI string length, then lexicographic order). The representative is always written in normalized (`/`-separated) form.

## Output

### BAM file (`--output`, required)

The input BAM with UMI tags updated. For each read whose UMI differs from its cluster's representative:
- The original UMI is moved to `--tag-orig` (default `OX`)
- The representative is written to `--tag-umi` (default `RX`)

Reads already matching the representative are unchanged. The output is coordinate-sorted via samtools sort.

### Molecule ID (`--write-mi`)

When `--write-mi` is set, each read receives an MI tag (name configurable via `--tag-mi`, default `MI`) with a two-level molecule group identifier: `mi_COMP.CLUST` (e.g., `mi_000001.001`). The first number identifies the read-overlap group (component) and the second identifies the UMI cluster within that component. All reads in the same UMI cluster share the same MI value. Clusters that were considered together during UMI clustering (i.e., in the same overlap group) share the same component number, making it possible to identify which clusters were compared against each other.

### UMI counts (`--summary-counts`, optional)

A tab-delimited BED6+ file with one row per UMI cluster per read-overlap-group. Coordinates are per-cluster (the bounding box of reads carrying that cluster's UMIs), not per-component. This file is not sorted by default. Columns:

| # | Column | Description |
|---|--------|-------------|
| 1 | chrom | Reference name |
| 2 | start | 0-based start of the cluster's reads |
| 3 | end | 0-based end (exclusive) of the cluster's reads |
| 4 | name | Molecule ID (e.g., `mi_000001.001` — component.cluster) |
| 5 | score | Total read count in the cluster |
| 6 | strand | `+`, `-`, or `.` |
| 7 | representative | Chosen representative UMI (slash-separated) |
| 8 | numUMIs | Number of distinct original UMIs in the cluster |
| 9 | maxEditDist | Max pairwise edit distance within the cluster (-1 if skipped for large clusters) |
| 10 | effectiveThreshold | Effective edit distance threshold after adaptive filtering (equals `--umi-edit-distance` when adaptive threshold is disabled or no distances were excluded) |
| 11 | umis | Comma-separated list of original UMIs |

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
| `--adaptive-alpha` | 0.10 | Maximum cumulative false positive rate per edit distance level |
| `--adaptive-threshold` | false | Discard edges at distances where cumulative random collisions exceed the FPR threshold |
| `--hp-dist` | false | HP-aware edit distance: one free HP indel per UMI segment (between separators) |
| `--ignore-refs` | | References to pass through without clustering (comma-separated) |
| `--junction-window` | 20 | Tolerance (bp) for matching junction positions and merging adjacent junctions |
| `--junction-match` | false | Require compatible splice junctions (CIGAR N ops) when grouping reads |
| `--match-one-end` | false | Match reads if EITHER 5' or 3' ends are within gap |
| `--no-strand` | false | Ignore strand when grouping reads |
| `--no-summary-counts-index` | false | Disable automatic tabix index generation for the summary counts file |
| `-o`, `--output` | (required) | Output BAM file path |
| `--overlap` | 50 | Maximum gap (bp) between reads to group them |
| `--region` | | Process only this region |
| `--write-mi` | false | Write the MI tag (molecule group ID) to output reads |
| `--tag-mi` | `MI` | SAM tag name for the molecule group ID written with `--write-mi` |
| `--tag-orig` | `OX` | SAM tag to store original UMI before correction |
| `--tag-umi` | `RX` | SAM tag containing UMI sequence |
| `-t`, `--threads` | 1 | Worker threads for UMI clustering |
| `--summary-counts` | | Write per-cluster UMI summary to this file |
| `--umi-cluster-method` | `adjacency` | Clustering method: `connected`,`adjacency`, `directional`, `tiered` |
| `--umi-edit-distance` | 3 | Maximum Levenshtein edit distance to cluster two UMIs |
`--whole-genome` | false | Process all UMIs as a single group (ignore coordinates) |

## References

- Smith T, Heger A, Sudbery I. "UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy." *Genome Research* 27(3):491-499, 2017. (Adjacency and directional clustering methods.)
