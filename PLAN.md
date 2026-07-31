# Roadmap

Forward-looking work that is **planned but not started**. Shipped behaviour is
documented in `README.md`, `CLAUDE.md` and `docs/`; this file exists so scoping work
isn't redone and so the reasoning behind deferred decisions survives.

---

## 1. Spanning VCF records are mis-indexed and mis-filtered (bug, blocks §2)

**Status: confirmed by reproduction, not yet fixed. Lives in `cghts/htsio/tabix`.**

A tabix region query misses any VCF record whose span extends past its first base.
Reproduced end to end with `cgkit`-built indexes, and isolated to the tabix layer
(`tabix.Reader.Query` alone, with no `varstore` involved):

| query | record it should find | found |
|---|---|---|
| `chr1:1400-1600` | `chr1:1000 N <DEL>` with `SVTYPE=DEL;END=2000` | **no** |
| `chr1:3100-3150` | `chr1:3000` 200 bp deletion, explicit `REF` | **no** |
| `chr1:4900-5100` | `chr1:5000 G>C` plain SNV (control) | yes |

### Where it actually goes wrong

**The reader is the primary cause.** `htsio/tabix/reader.go:502-509` computes
`end := beg + 1` for VCF (overridden only `if meta.ColEnd != 0`, and the VCF preset
sets `colEnd = 0`), and `reader.go:458` filters on it: `if rec.End <= start`. A record
whose true span reaches into the query region is therefore discarded by the overlap
test regardless of what the index says. This is why the reproduction above fails, and
it is **independent of how the file was indexed** — a bcftools-built `.tbi` would not
help.

**The writer is a second, distance-dependent cause.** `writer.go:513-520` bins by
`Reg2Bin(l.start, recEnd)` with the same `recEnd = l.start + 1`. But tabix's smallest
bin covers 16 kb, so for the cases above a point and a span land in the *same* bin
(both `4681`) and the chunk is still examined:

| interval | bin |
|---|---|
| 200 bp del @3000 binned as a point | 4681 |
| same, binned by `len(REF)` | 4681 |
| query `3100-3150` | 4681 |
| a 100 kb del @1000 binned as a point | 4681 |
| query `90000-90100` | **4686** |

So point-binning only loses a record once its span crosses a 16 kb window — a large
SV or a long gVCF block. Fixing the reader alone would repair the reproduced cases;
fixing both is required for spans over 16 kb.

### The span's two sources

- **`len(REF)`.** htslib indexes a VCF record's end as `beg + len(REF)`; cghts uses
  `beg + 1`. This affects **plain long deletions and MNPs**, with no symbolic ALT and
  no `INFO/END` anywhere — the `chr1:3100-3150` reproduction is exactly this case.
- **`INFO/END`**, for symbolic ALTs (`<DEL>`, `<DUP>`, `<CNV>`) and gVCF reference
  blocks.

> **Open question — worth settling before implementing.** Marcus's read is that tabix
> cannot carry an end differing from a column value, which would mean the `INFO/END`
> half is not htslib-compatible behaviour to match. My reading of htslib's `tbx.c`
> `get_intv` is that the `.tbi` `format` field selects preset-specific end derivation:
> `TBX_SAM` computes the end from the **CIGAR** (not a column under any reading), and
> `TBX_VCF` uses `len(REF)` and scans INFO for `END=`. `col_end = 0` for VCF is what
> makes that logic necessary rather than what forbids it. **Unverified here** — no
> `tabix`/`bgzip`/`bcftools` in this environment, so this rests on source knowledge,
> not a measurement. Settle it by indexing the same SV VCF with htslib and diffing
> query results. Note the `len(REF)` half stands either way, so §1 is a real bug
> regardless of how this resolves.

This is not only a gVCF prerequisite. It is a present-day correctness bug for
long-indel VCFs — and, if the `INFO/END` reading holds, for structural-variant VCFs,
which are a supported input class here (`vcf-svtofasta`, `vcf-tobedpe`). Every path
writing a VCF `.tbi` is affected: `vcf-filter --tbi`, `vcf-varquery --tbi`,
`vcf-split`, `tab-index`, `tab-sort`.

Worth noting what this is *not*: the end-of-record machinery already exists and the
BED preset exercises it. The fix supplies a VCF-shaped source for the value rather
than building interval support from scratch.

**Fix order:** reader first — it is the cause that always bites, and it is what makes
externally-built indexes usable. Then the writer, for spans over 16 kb. Then
differential-test against an htslib-produced `.tbi` on the same file. Note
`htsio/tabix/` has had other work in flight — check for overlap before starting.

---

## 2. gVCF support

**Status: planned, not started.** Nothing in cgkit or cghts references gVCF today.
Blocked on §1.

### Why it's worth doing

A plain VCF asserts nothing about the positions between its records. That single fact
constrains the entire genotype store: the sites catalog is the exact boundary of what
is answerable, an off-catalog locus is `not_assayed` for every sample, and
`non_carrier` is *inferred* from adjacent variant sites rather than observed.

A gVCF differs in kind. Its reference blocks (`<NON_REF>`/`.` ALT, `INFO/END`,
`MIN_DP`, `RGQ`) are **positive statements about spans**. Three things follow that
nothing else provides:

1. **Observed non-carriers.** A reference block covering the position with `MIN_DP`
   above threshold is direct evidence, not an inference from neighbours.
2. **Off-catalog answers.** `varstore.Classify` returns early on the catalog check
   specifically so a sites-store cannot overreach; that check is already written to
   admit a blocks-store.
3. **A query-time depth threshold.** Today `--min-dp` is baked into the callable runs
   at conversion, so a query using a different value silently disagrees with the
   store. A gVCF carries banded `MIN_DP` per block, so the number can come from the
   caller instead.

The store format already anticipates this: `varstore.SpansBlocks` exists as a
`SpanSemantics` value and is documented as not yet produced. Anything without the
metadata key reads as `SpansSites`, the conservative interpretation, so adding blocks
reinterprets nothing already written.

### Sequencing

1. Fix §1 in cghts (tabix writer + reader).
2. cghts `varstore`: produce `SpansBlocks` stores; make `Classify` and reference-call
   reconstruction treat block coverage as observation. The query paths already branch
   on `SpanSemantics`.
3. cgkit: `vcf-toparquet` accepting gVCF input; a query-time `--min-dp` no longer
   required to match the conversion value.
4. Release cghts, then `make bump-cghts` here (see `CLAUDE.md`).

### Traps already identified

- **The cross-backend equivalence tests will not cover this.** A blocks-store
  answering off-catalog positions has no VCF-backed equivalent to compare against, so
  the mechanism that caught the sibling-allele bug is blind here. It needs independent
  ground truth (bcftools/GATK on the same gVCF).
- **Run intervals in a `SpansSites` store must keep meaning what they mean.** The risk
  is a change that makes `regions.parquet` look like coverage for *all* stores,
  retroactively licensing claims the source never made.
- **Do not "upgrade" existing stores by inference.** `SpansSites` stays the default
  for anything lacking the metadata key.

---

## 3. Bytes-read / row-groups-decoded counters (optional)

`internal/cmd/vcfcmd/bench_test.go` measures wall time, which is hardware-bound — the
durable claim it supports is the *ratio* (bulk flat vs per-locus linear), not the
milliseconds. Counters inside `varstore` would give numbers that survive a machine
change. Not required by anything.

---

## Decided against — do not re-litigate

- **`vcf-gtmatrix` as a separate command.** A VCF with GT already *is* the wide
  matrix; `vcf-varquery --format vcf` covers it.
- **`--classify`, `--min-gq`, `--region`.** The genotype in the output already says
  which case occurred, so there is no state column to interpret. `--min-gq` was
  *broken*, not merely redundant: callable runs are built from depth alone, so a store
  retains no GQ for a genotype it never wrote down and could not gate a reference call
  where a VCF would — the backends silently disagreed. `--region` folded into
  `--variant`.
- **Page-level pruning in `varstore`.** Retired by measurement: the bulk path is flat
  out to 1000 targets, so page pruning would only shave the single-locus case — not
  worth added risk in code that must never skip a matching row group.
- **Multi-locus tabix seeking in `VcfStore`.** Only pays off for panel queries against
  plain VCFs rather than converting first, which is not the workflow the store serves.
- **A floating query-time `--min-dp` without gVCF.** Per-run `MinDP` alone is
  pathological (one DP-5 site poisons a thousand DP-30 sites in the same run); DP
  banding fragments toward one row per site; per-site DP is just the dense matrix
  again. The number has to come from the source — hence §2.
- **A sample-sorted copy of `calls.parquet`.** Decode dominates, not layout, and it
  would break the `(chrom, pos)` sort that locus pruning depends on.
