# Roadmap

Forward-looking work that is **planned but not started**. Shipped behaviour is
documented in `README.md`, `CLAUDE.md` and `docs/`; this file exists so scoping work
isn't redone and so the reasoning behind deferred decisions survives.

---

## 1. cghts's tabix layer treats every VCF record as one base wide (bug, blocks §2)

**Status: confirmed against htslib source. Lives in `cghts/htsio/tabix`.**

A region query misses any VCF record whose span extends past its first base.
Reproduced end to end with cgkit-built indexes, then isolated to `tabix.Reader.Query`
with no `varstore` involved:

| query | record it should find | found |
|---|---|---|
| `chr1:3100-3150` | `chr1:3000`, 200 bp deletion, explicit `REF` | **no** |
| `chr1:1400-1600` | `chr1:1000 N <DEL>`, `SVTYPE=DEL;END=2000` | **no** |
| `chr1:4900-5100` | `chr1:5000 G>C` plain SNV (control) | yes |

### What the format actually specifies

The `.tbi` header's `col_end` is `0` for the VCF preset, which reads naturally as "VCF
records have no end". The TBI spec (`hts-specs/tabix.tex`) has two separate rules, and
the distinction is the whole thing:

> For the SAM format, the end of a region equals `POS` plus the reference length in the
> alignment, inferred from `CIGAR`. For the VCF format, the end of a region equals `POS`
> plus the size of the deletion.
>
> Field `col_beg` may equal `col_end`, and in this case, the end of a region is
> `end`=`beg+1`.

The second rule is conditioned on **`col_beg == col_end`**, not on `col_end == 0`. For
VCF, `tbx_conf_vcf = { TBX_VCF, 1, 2, 0, '#', 0 }` gives `col_beg = 2`, `col_end = 0`,
so they are unequal and that rule does not apply; the first rule does. htslib matches:
`tbx.c:136` reads `if (conf->bc <= conf->ec) intv->end = intv->beg;`, which for VCF is
`2 <= 0` — false — so the end is deliberately left unset at the POS column and filled
in from REF later.

Both bullets are byte-identical to the original 2010 commit that added the spec
(`385b4743`); the file has had five commits ever, none touching this text. So this is
not a rule that changed under anyone.

**Two fair criticisms of the spec, since they explain why this is easy to get wrong.**
It states a format-specific coordinate rule as an unnumbered note rather than in the
`col_end` field description, where a reader looking up `col_end = 0` would actually
find it. And "the size of the deletion" is vague — it means the reference span,
`len(REF)`, which also covers MNPs and is `1` for an insertion.

**One caveat worth keeping straight.** The spec documents only the `len(REF)` case. It
says nothing about `INFO/END`, `SVLEN` or `FORMAT/LEN`; those are htslib behaviour
beyond the written spec. So matching the spec fixes long indels, while matching htslib
is what makes SV and gVCF work — and only the former can be justified from the spec
alone. Anything reading a `.tbi` written by htslib has to match the implementation, not
the document.

### What htslib implements

`col_end` is consulted only for `TBX_GENERIC`; each other preset derives its own end.
From `tbx.c` on `develop`, read directly:

| preset | line | end is computed from |
|---|---|---|
| `TBX_GENERIC` | 153-156 | the `col_end` **column** — the only preset that reads one |
| `TBX_SAM` | 170 | `intv->end = intv->beg + l`, `l` summed over **CIGAR** `M`/`D`/`N` |
| `TBX_VCF` | 174 | `intv->end = intv->beg + (i - b)` — i.e. `len(REF)` |
| `TBX_VCF` | 230 | `intv->end = end` from **`INFO/END=`** (ignored, with a warning, if `END <= POS`) |
| `TBX_VCF` | 206-300 | `INFO/SVLEN` for symbolic ALTs, and `FORMAT/LEN` for `<*>`/`<NON_REF>` — the ALT scan at line 197 is commented `//note gvcf` |
| `TBX_VCF` | 304-309 | reconciles them: `max(reflen, svlen, fmtlen) + beg`, then keeps whichever is larger, that or the `INFO/END` value |

So a calculated end is the rule for the structured presets, not the exception —
`TBX_GENERIC` is the only one that reads an end column at all. SAM is the clearest
case: a CIGAR-derived end is not a column value under any reading. Line 310 even notes
the calculation is kept in sync with `vcf.c:get_rlen`.

**This applies to `.tbi` specifically, not just `.csi`.** `tbx_index()` line 493 takes
`min_shift == 0` → `n_lvls = 5, fmt = HTS_FMT_TBI`; line 532 then calls
`hts_idx_push(tbx->idx, intv.tid, intv.beg, intv.end, …)` with that same computed
`intv.end` in **both** the TBI and CSI branches. The choice only changes bin geometry
(`min_shift`, `n_lvls`) and the on-disk container — never the interval being indexed.
`tbx_readrec` (line 371) hands the same `intv.end` to the iterator, so the read-side
overlap test uses it too.

### Where cghts diverges

cghts implemented the `TBX_GENERIC` column path and never added the preset-specific
derivation, in either direction:

- **Reader — the cause that always bites.** `reader.go:502-509` computes
  `end := beg + 1`, overridden only `if meta.ColEnd != 0`, and `reader.go:458` filters
  on it (`if rec.End <= start`). So a spanning record is dropped by the overlap test
  **regardless of how the file was indexed** — an htslib-built `.tbi` would not help.
- **Writer — secondary, distance-dependent.** `writer.go:513-520` bins by
  `Reg2Bin(l.start, recEnd)` with the same understated `recEnd`. But TBI's smallest bin
  covers 16 kb, so a point and a short span usually land in the *same* bin and the
  chunk is still examined:

  | interval | bin |
  |---|---|
  | 200 bp del @3000 as a point | 4681 |
  | same, by `len(REF)` | 4681 |
  | query `3100-3150` | 4681 |
  | 100 kb del @1000 as a point | 4681 |
  | query `90000-90100` | **4686** |

  So point-binning only loses a record once its span crosses a 16 kb window — a large
  SV, or a long gVCF block.

This is not only a gVCF prerequisite. It is a present-day correctness bug for
long-indel and structural-variant VCFs, which are a supported input class here
(`vcf-svtofasta`, `vcf-tobedpe`). Every path writing or reading a VCF `.tbi` is
affected: `vcf-filter --tbi`, `vcf-varquery --tbi`, `vcf-split`, `tab-index`,
`tab-sort`.

Worth noting what this is *not*: the end-of-record machinery already exists and the BED
preset exercises it. The fix supplies a VCF-shaped source for the value rather than
building interval support from scratch.

**Fix order:** reader first — it is the cause that always bites, and it is what makes
externally-built indexes usable. Then the writer, for spans over 16 kb. Match htslib's
precedence (`len(REF)`, then `INFO/END` when greater than POS, then `SVLEN`/`LEN`), and
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
