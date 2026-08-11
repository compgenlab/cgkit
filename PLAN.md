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

### Conversion is a cohort operation — deferred to the next cgkit release

Not a stopgap, and not a single-file feature we simply have not written yet. **The
unit of gVCF conversion is a *set* of gVCFs, converted together.** A gVCF is
single-sample by construction, so converting one produces a one-sample store — which
answers nothing a direct `vcf-varquery` against that gVCF does not already answer,
faster and with no conversion step. The value of a store is the cohort: many samples
at one locus, which is exactly what one gVCF cannot give.

That reframes the work. It is not "teach the converter about blocks" plus "later,
handle several inputs" — handling several inputs *is* the feature, and the block
handling is a detail inside it.

What it needs, in rough order of difficulty:

- **An N-way coordinate merge.** N single-sample gVCFs all start at chr1:1, so
  today's `checkOrder` rejects them outright, and simply concatenating would break
  the `(chrom, pos)` sort that every row-group pruning bound depends on.
- **A union roster.** The converter requires every input to carry the *same* samples;
  for gVCFs each input brings a different one. The mechanical part is small (both
  positional uses already route through `perm`), and the new error is the opposite of
  today's: the same sample appearing in two inputs, which would double-count.
- **A decision on AC/AN.** A site row is written the moment its record is read, and
  there is no cross-input accumulator. Merging allele counts across N gVCFs at a
  shared locus is joint genotyping. Either compute them over the merged cohort or
  refuse to emit them — but not the current behaviour, which would report a
  per-sample count as if it were a cohort count.
- **The block schema**: block spans in `regions.parquet` with a per-block `MIN_DP`
  column, `SpansBlocks` set, `Classify`'s off-catalog branch wired to the regions
  scan (it already tests `s.spans != SpansBlocks`), and `callsWithRef` given a second
  emission source — it is driven by the sites catalog, so it emits nothing
  off-catalog no matter what the metadata says.

Until then `vcf-toparquet` refuses a gVCF, which is the right answer rather than a
placeholder: the store it would write is wrong in ways nothing downstream can detect.

The refusal keys on **a reference block record**, not on the header. It once used
`isGvcfHeader` and rejected ordinary cohort VCFs for it: a DRAGEN msVCF is a
joint-genotyped callset — precisely this command's input — and still carries the
`##ALT=<ID=NON_REF>` line it inherited from the gVCFs it was built from, as do many
GATK-derived callsets. `##ALT=<ID=*>` is worse than incidental: VCF 4.5 gives the
gVCF unspecified allele and the ordinary spanning-deletion allele the same ALT ID,
so no header test can separate them even in principle. `vcf-strip` still has to
decide from the header, because it writes its output header before reading a record;
`vcf-toparquet` does not, and a failed conversion `Discard()`s, so reading first
costs nothing recoverable. A mixed `G,<NON_REF>` record is a variant record and
converts, with the block allele masked out the way any non-focal alternate is.

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

## 4. Indel normalization — not started, and mostly not ours

`vcf-toparquet` splits multiallelics but does **not** normalize allele
representation. Splitting and normalization are routinely conflated because one
`bcftools norm` invocation does both; they are separate operations and only the
first is handled here. `bcftools norm -m -both` — a common pre-step — does the
half cgkit already does and none of the half it does not.

The exposure is exact-locus matching. `SameLocus` canonicalizes the chromosome
name and then compares `Pos`/`Ref`/`Alt` byte-for-byte
(`cghts/varstore/store.go:426`), so an indel stored as the source wrote it and
queried in a different-but-equivalent representation is a miss — no rows, and the
target-not-matched warning is the only signal. Span queries (`chr1:1000-2000`, a
whole contig, BED) are immune: overlap comes from `pos`/`ref_end` and never touches
the alleles. SNVs are immune in practice. The exposure is indels reached by exact
locus, and any join between a store and an externally-normalized site list.

Normalization has two halves with very different costs.

**Trimming is reference-free.** Right-trim while all alleles have length ≥ 2 and
share a last base; left-trim while all alleles have length ≥ 2 and share a first
base, incrementing `POS`. So `chr1:100 CTT>CT` → `chr1:100 CT>C`, and
`chr1:100 GAT>GA` → `chr1:101 AT>A`.

**Left-shifting needs the FASTA**, because the shift is entirely a property of the
flanking reference. With reference `97:C 98:A 99:A 100:A 101:A 102:A 103:G`, a
one-base deletion written `chr1:101 AA>A` left-aligns to `chr1:97 CA>C`. Nothing in
the record predicts that distance: it is the length of the perfect tandem repeat
containing the indel — tens of bases for a homopolymer, hundreds for an STR,
kilobases at an expanded repeat locus.

### Three things it breaks

1. **Sort order.** `checkOrder` errors on `pos < lastPos`
   (`internal/cmd/vcfcmd/vcf_toparquet.go:833`), which is the correct failure but
   means normalization cannot be a filter inside the record loop. It needs a
   reorder buffer, and the window cannot be computed — only chosen, with an
   explicit error when a record shifts past it. `bcftools norm -w/--site-win`
   (default 1000) is the precedent worth copying rather than reinventing.
   Trimming needs one too, just a small one: left-trim moves `POS` *forward* by up
   to `len(REF)-1`, so record N can trim past record N+1.
2. **Row-group pruning.** Tight per-group min/max on `pos` is what gives the
   measured 3.05x, so reordering has to happen before the writer sees a row.
3. **`RecordKey`.** Sibling exclusion for `--hom-ref` keys on `{Chrom, Pos, Ref}`
   derived from the locus (`cghts/varstore/store.go:52`): a `0/2` sample is kept out
   of the allele-1 row because both split loci share a key. Normalizing *after* the
   split breaks that, under trimming alone — `chr1:100 GA → GAA, G` splits to
   `chr1:100 G>GA` and `chr1:100 GA>G`, same position, different `REF`, different
   key — and a `0/2` sample silently reappears as `0/0`. Preserving it means
   carrying the source record identity explicitly, which is a schema change to
   `calls`/`sites`, not a flag. **This is the real cost of the feature**, and it is
   easy to miss when scoping it as "call a normalize function per record".

### The hazard that argues for doing something now

A missing FASTA fails loudly. A FASTA for the **wrong build** does not: shifting
against the wrong bases yields well-formed records at confidently wrong positions,
and the store then answers every query in coordinates that mean something else.
Same failure class `--meta-reference` exists for, but worse — metadata records a
claim, this manufactures data.

The guard is cheap and header-only, so it should be mandatory rather than opt-in:
compare `(*vcf.VcfHeader).ContigDef(id).Length` (`cghts/vcf/header.go:178`, `-1`
when absent) against `seqio.ReferenceReader.SequenceLength(name)`
(`cghts/seqio/reference.go:24`), matching names through `CanonicalContig` since a
`chr1` FASTA against a `1` VCF is routine. No records read, so it costs nothing
against a 200 GB input and it catches GRCh37-vs-38 immediately.

### If any of this gets built, build this first

A **reference-free detector**: per split allele, flag `len(ref) ≥ 2 && len(alt) ≥ 2
&& (ref[0] == alt[0] || ref[len-1] == alt[len-1])`, count the hits, and report them
in the `-v` conversion summary beside the DP/GQ/AD census, plus a manifest field.
It flags `GAT>GA` and `CTT>CT` and leaves `GA>G` and `A>T` alone. It must run
**after** the split, not before: `GA → GAA, G` is minimal as a multiallelic (the `G`
allele is length 1, so no joint trim applies) while its split components are not.

That is nearly free, needs no FASTA, no flag, no reorder buffer and no schema
change, and it converts a silent mismatch into a stated fact about the store. Full
left-alignment is the project; detection is the part that pays for itself.

### The counter-argument, recorded so it is not rediscovered

`bcftools norm -f ref.fa` already does this correctly, is universally available, and
runs upstream of conversion. The case for owning it in cgkit is **not** capability —
it is that a store should be self-describing about whether it was normalized and
against what, so two stores can be known to be joinable. That is a metadata goal
that the detector above serves at a fraction of the cost. Do not build the
realignment for the capability alone.

Also note that pre-normalizing with `-m -both` is now actively counterproductive:
cgkit tallies `AC`/`AN` from the record's raw GTs before recoding
(`AddAlleleCounts`, `cghts/varstore/vcfrecord.go:139`), so it sees both alleles of a
`1/2`. A record arriving already split with `.` padding has lost that, and `AN` —
and therefore `AF` in `--format vcf` — comes out low at every multiallelic site.
Recommend `bcftools norm -f ref.fa` with no `-m`.

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
