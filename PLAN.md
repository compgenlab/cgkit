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

## 2. gVCF support — NOT started; §1 was the prerequisite, not the feature

**Status: the plumbing is in place, gVCF as a data source is not.** §1 means a gVCF
can now be indexed, region-queried and passed through `vcf-strip` without being
silently destroyed. It does not mean the store understands reference blocks.

Converting a gVCF today **appears to work and is wrong**. From
`testdata/sample.g.vcf`:

```
$ vcf-toparquet --out g/ sample.g.vcf     # 2 blocks + 1 variant
wrote g/: 1 calls, 4 sites, 1 callable runs over 1 samples

$ vcf-varquery --variant chr1 --hom-ref g/
chr1  100   A  <NON_REF>  S1  0/0  .   10  ...
chr1  5001  A  G          S1  0/1  44  44  ...
chr1  5002  C  <NON_REF>  S1  0/0  .   10  ...
```

Four things wrong there, and none is a crash:

1. **`<NON_REF>` enters the sites catalog as a variant.** A reference block becomes
   a pseudo-site with an ALT of `<NON_REF>`, so the catalog — which is supposed to
   be the exact boundary of what is answerable — is polluted with non-variants, and
   `AC`/`AN` count a block allele as an allele.
2. **`span_semantics` is still `sites`.** `varstore.SpansBlocks` exists and
   `Classify` already branches on it, but nothing sets it, so an off-catalog
   position is `not_assayed` even though the gVCF explicitly asserted coverage
   there.
3. **`MIN_DP` is ignored.** `min_dp` reports the conversion `--min-dp` (10), not the
   block's own `MIN_DP` (28). This is what a query-time depth threshold needs.
4. **The block reads as a variant call.** Because §1 landed, a query *inside* a
   block does now return it by overlap — `--variant chr1:2000-2100` yields the row
   at 100 — which is accidentally useful and semantically wrong: it is reported as
   a `<NON_REF>` variant rather than as reference across a span.

So the work remaining is about **meaning**, not retrieval:

- `vcf-toparquet` recognises a gVCF (`isGvcfHeader` already exists in `vcfcmd`, used
  only by `vcf-strip`) and treats a reference block as coverage rather than as a
  variant: no catalog site, no allele counts, no ALT call.
- `regions.parquet` records block spans, and the store declares `SpansBlocks`.
- `Classify` and reference-call reconstruction treat block coverage as an
  observation, which is the only way `non_carrier` becomes observed rather than
  inferred from adjacent variant sites.
- `MIN_DP`/`RGQ` drive the gate, enabling a query-time `--min-dp` that need not
  match the conversion value.
- A mixed record (`G,<NON_REF>`) keeps its real allele and drops the block one.
  `vcf.IsRefBlockAlt` exists for exactly this.

### Why it is worth doing

A plain VCF asserts nothing about the positions between its records, so a missing
row cannot be told from an unsequenced one. A gVCF's blocks are positive statements
about spans, which is what turns "absence of evidence" into "evidence of absence" —
the difference between *"not reported"* and *"confidently reference at depth 30"*.
That distinction is the whole question for a polygenic score, where treating an
absent variant as `0/0` biases the result and treating it as missing drops a term,
with nothing reporting either.

### Traps already identified

- **The cross-backend equivalence tests cannot cover this.** A blocks-store
  answering off-catalog has no VCF-backed equivalent to compare against, so the
  mechanism that caught the sibling-allele bug is blind. It needs independent ground
  truth (bcftools/GATK on the same gVCF).
- **Run intervals in a `SpansSites` store must keep meaning what they mean.** The
  risk is a change that makes `regions.parquet` look like coverage for *all* stores,
  retroactively licensing claims the source never made.
- **`SpansSites` stays the default** for anything lacking the key. Do not upgrade
  old stores by inference.
- Conversion reads sequentially, so it never needed §1. A `--region`-bounded
  conversion would, and cannot simply seek to the region: a block overlapping the
  region's start begins before it.

### Not in scope

Writing gVCF output; emitting reference blocks from `vcf-tobed` or
`vcf-varquery --format vcf`; joint genotyping or merging N gVCFs (that is
GATK/GenomicsDB's job — we consume its output); per-base depth.

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
