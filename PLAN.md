# Roadmap

Forward-looking work that is **planned but not started**. Shipped behaviour is
documented in `README.md`, `CLAUDE.md` and `docs/`; this file exists so scoping work
isn't redone and so the reasoning behind deferred decisions survives.

---

## 1. Record spans — DONE in cghts v0.7.7

Kept as a record of what the fix was, since the reasoning is load-bearing for §2.

A tabix region query used to miss any VCF record reaching past its first base: a
long deletion, a symbolic `<DEL>`, a gVCF reference block. `chr1:3100-3150` over a
200 bp deletion at `chr1:3000` returned nothing.

The TBI spec has two coordinate rules and VCF falls under the first — "for the VCF
format, the end of a region equals `POS` plus the size of the deletion" — while
`end = beg+1` applies only when `col_beg == col_end`. `col_end = 0` is not that
case. htslib agrees in `tbx_parse1`: `TBX_GENERIC` is the only preset that reads an
end column; `TBX_SAM` sums CIGAR `M`/`D`/`N`, and `TBX_VCF` uses `len(REF)`,
`INFO/END`, `INFO/SVLEN` and `FORMAT/LEN`. Both bullets are byte-identical to the
2010 commit that added the spec, so this was never a rule that changed.

What shipped:

- `internal/vcfspan` holds the precedence once, shared by `htsio/tabix` (split
  lines) and `(*vcf.VcfRecord).RefSpan` (parsed records), because tabix cannot
  import `vcf`. A drift test holds both to one corpus.
- The `.tbi` header now records its preset. It was hardcoded to generic, so nothing
  could apply VCF rules on read. **Earlier indexes still read as generic and must
  be rewritten.**
- `varstore` selects by overlap via a `ref_end` column; pruning's lower bound comes
  from `max(ref_end)`, since a matching record may start below a row group's
  minimum.
- `vcf-strip` refuses to drop `INFO/END` from a gVCF without `--force-end`;
  `vcf-tobed` reports true spans and skips reference blocks; `vcf-toparquet`
  records the span.

## 2. gVCF support — querying DONE; conversion deliberately refused

**`vcf-varquery` reads gVCFs.** A reference block asserts that a sample was called
reference across a span at at least `MIN_DP`, and that is now what it reports:

```
$ vcf-varquery --variant chr1:2000-2100 --hom-ref --min-dp 10 sample.g.vcf.gz
chrom  pos  ref  alt  sample  gt   dp  min_dp  ad_ref  ad_alt  gq
chr1   100  A    .    S1      0/0  .   28      .       .       60
```

`alt` is `.` because a block names no alternate; `min_dp` is the block's own
`MIN_DP`, not the query threshold; `dp` is `.` because a block measures no single
base; `gq` carries `RGQ`. One row per block per sample — never one per base.

The boundary still holds. A **gap** between blocks was never reported, so it stays
unanswerable, and a block below the gate is excluded. `Gate.Admits` prefers `MinDP`
over `DP` for that reason: `Missing` passes every gate by design, so a block-derived
call with no `DP` would otherwise be admitted unconditionally — the gate silently
doing nothing on exactly the state gVCF makes trustworthy.

**`vcf-toparquet` refuses a gVCF.** It used to convert one happily and produce a
store wrong in three ways nothing downstream could detect: `<NON_REF>` in the sites
catalog, `AC`/`AN` counting a block allele as an allele, and every block reduced to
its first base.

### What conversion would still need

Deferred, not forgotten. A blocks store needs block spans in `regions.parquet` with
a per-block `MIN_DP` column, `SpansBlocks` set, `Classify`'s off-catalog branch wired
to the regions scan (it already tests `s.spans != SpansBlocks`), and `callsWithRef`
given a second emission source — it is driven by the sites catalog, so it emits
nothing off-catalog no matter what the metadata says.

Merging N gVCFs is a separate problem again: they are single-sample, so a cohort is N
inputs each bringing a *different* sample. That inverts the converter's same-samples
rule, breaks `checkOrder` (they all start at chr1:1), needs an N-way merge to keep
the `(chrom,pos)` sort every pruning bound depends on, and raises AC/AN denominators
across inputs — which is joint genotyping, and out of scope.

### The testing note that changed

An earlier version of this plan recorded that a blocks store answering off-catalog
would have no VCF-backed equivalent, leaving the cross-backend harness structurally
blind. Making `VcfStore` the *first* block-aware backend inverts that: a gVCF read
directly is now the ground truth a future blocks store gets checked against.

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
