# vcf-toparquet: the Parquet genotype store

A columnar, sparse on-disk format for cohort genotypes, written by `cgkit vcf-toparquet` and queried by `cgkit vcf-varquery`. It answers "who carries this variant" and "what does this subject carry" without reading every genotype in the callset, and — unlike a plain VCF — it keeps enough context to tell a confidently-called reference genotype apart from a position that was never assayed.

## Usage

```
cgkit vcf-toparquet --out DIR <input.vcf> [input2.vcf ...]
cgkit vcf-varquery [--variant LOCUS | --sample NAME] <input.vcf | store>
cgkit vcf-varsummary <input.vcf | store>
```

## What it is

**A store is a directory**, and its members are one inseparable set — the reasoning is in [Why three files](#why-three-files-and-not-one).

```
cohort/
  calls.parquet      one row per ALT-carrying genotype
  sites.parquet      one row per interrogated site, with allele counts
  regions.parquet    one row per run of called sites, per sample
  manifest.json.gz   written last; the store is unreadable without it
```

A trailing `/` on `--out` is optional and means nothing: `cohort` and `cohort/`
name the same store, and `vcf-varquery` also accepts any member path within it.
There used to be a second, filename-prefix form (`cohort.calls.parquet`); it is
gone, because keeping both meant every path decision first had to guess which
was meant — and the guess required a filesystem check, which cannot work for a
remote locator at all.

A directory is also what every Parquet tool expects, so DuckDB or pyarrow can be
pointed straight at `cohort/calls.parquet`.

> **Migrating.** Two breaking changes land together: stores are directories, and
> the manifest is required. A store written by an earlier cgkit has neither and
> must be re-converted. `cgkit vcf-varsummary <store>` reports why one is
> refused.

### Metadata: what the store *is*

Everything above records how a store was **made** — the command, the inputs, the
sample roster, the moment it ran. None of it says what the store **holds**.

A whole-genome callset ships as one VCF per chromosome and converts into a single
store. That store can name its 24 input filenames, but not the release they
collectively are — and those filenames stop identifying anything the moment the
store is copied to another machine. The assembly is the same gap in a more
dangerous form: the `##contig` lines carry lengths that *imply* GRCh37 or
GRCh38, but implying is not declaring, and a store read against the wrong
assumption does not fail. It answers, with coordinates that mean something else.

```
cgkit vcf-toparquet \
  --meta-dataset 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05 \
  --meta-reference GRCh38 \
  --meta-caller 'GATK 4.2.6.1' \
  --meta cohort=phase3 \
  --out cohort chr*.vcf.gz
```

Seven keys have an agreed meaning — `dataset`, `reference`, `caller`,
`accession`, `url`, `version`, `description` — each with a `--meta-<key>` flag.
`--meta KEY=VALUE` records anything else and repeats.

It lands in two places, the same two that already carry provenance: the
manifest, as a `meta` object, and the calls file's Parquet key/value metadata,
namespaced under `cgkit.meta.` so a supplied key can never shadow one the writer
vouches for (`--meta source=...` does not become `cgkit.source`).

```jsonc
// manifest.json.gz
{
  "sources": ["chr1.vcf.gz", "..."],   // how it was made
  "meta": {                            // what it is
    "dataset":   "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05",
    "reference": "GRCh38",
    "caller":    "GATK 4.2.6.1",
    "cohort":    "phase3"
  }
}
```

Values are recorded **verbatim and never validated** — cgkit cannot know whether
`GRCh38` is true, and normalizing it would quietly turn your claim into cgkit's.
Keys *are* validated, to lowercase `[a-z0-9_-]`, because a key has to survive a
Parquet metadata key, a JSON member and a `grep`.

A key given twice, by either spelling, is an **error** rather than last-wins:
which of two conflicting claims about a store gets recorded should not depend on
the order the flags were typed.

Metadata is optional, and a store converted without any omits the field
entirely rather than writing `{}` — absent means "not stated", never "stated as
nothing". `vcf-varsummary` prints all of it, `--format json` emits it verbatim,
and `vcf-varquery -v` reports `dataset` and `reference`, the two that change how
a result should be read.

### Schemas

`calls.parquet` — one row per (sample, site) where the sample carries the ALT allele. Records are normalized to one variant per row, so `alt` is always a single allele.

| column | type | notes |
|---|---|---|
| `sample_id` | string (dict) | bloom-filtered |
| `chrom` | string (dict) | the source's own spelling, preserved |
| `pos` | int32 | 1-based |
| `ref`, `alt` | string (dict) | one ALT allele per row |
| `gt` | string (dict) | recoded biallelic genotype (see [normalization](#normalization)) |
| `dp` | int32 | `-1` when absent |
| `ad_ref`, `ad_alt` | int32 | per-allele, not summed; `-1` when absent |
| `gq` | int32 | `-1` when absent |

`sites.parquet` — one row per interrogated site, independent of any sample. A site with zero carriers still gets a row; that is the entire point of the file.

| column | type | notes |
|---|---|---|
| `chrom`, `pos`, `ref`, `alt` | | the site's identity |
| `ac` | int32 | ALT **alleles** observed |
| `an` | int32 | called **alleles** at the site |
| `n_carriers` | int32 | **samples** with ≥1 copy of this ALT |
| `n_called` | int32 | **samples** called and passing `--min-dp` |
| `n_lowdp` | int32 | **samples** failing `--min-dp` here |

`regions.parquet` — runs of consecutive catalog sites at which one sample was called at adequate depth.

| column | type | notes |
|---|---|---|
| `sample_id` | string (dict) | bloom-filtered |
| `chrom` | string (dict) | runs never span chromosomes |
| `start`, `end` | int32 | first and last **catalog site position** in the run |
| `n_sites` | int32 | how many catalog sites the run covers |

**`start` and `end` are site positions, not an interval of covered bases.** See [runs are not coverage](#runs-are-not-coverage).

Absent integer fields are stored in-band as `-1`, exposed as `varstore.Missing`. VCFs routinely omit DP/GQ/AD — a GT-only phased reference panel has none of them — and the columns are kept non-optional so reads stay a flat scan. Callers must test against `Missing` before using a value; a naive comparison would read it as an extremely low quality score.

## Why it exists

A VCF is row-major over samples: one line per site, one column per sample. That layout has a cost that grows with the *cohort*, not with the question. Asking "does HG00096 carry anything in BRCA1" means parsing every genotype of every other sample at every site in the range, because they are interleaved on the same lines. Asking "who carries chr17:43045642 G>C" in a 3,202-sample callset means decompressing a line that is mostly the 3,201 answers you did not ask for.

Two properties fix that:

**Columnar.** Parquet stores each column contiguously, so a query touches only the columns it names, and per-column-chunk statistics let a reader skip whole blocks. `vcf-toparquet` writes rows in coordinate order and declares that ordering, which makes the per-row-group `min`/`max` on `pos` tight and non-overlapping — so a locus lookup can skip every block that cannot contain it. See [row-group pruning](#row-group-pruning).

**Sparse.** Only ALT-carrying genotypes are stored. Most genotypes in a cohort callset are homozygous reference, and writing them down costs a great deal to say very little. On the 1000 Genomes 30x BRCA1 slice — 3,747 sites × 3,202 samples ≈ 12.0M genotypes — the store holds **439,183 calls, about 3.7%**.

Sparsity is what creates the format's one real problem, and the next section is about solving it.

## Why three files, and not one

Drop the other two files and the store can no longer distinguish these two situations:

- a sample was sequenced, called at good depth, and is **reference** at this position
- a sample was **never assayed** here at all

Both look identical: no row in `calls.parquet`. Collapsing them into "not a carrier" silently inflates the denominator of every cohort query, every allele frequency, and every association test. So the store has to be able to tell four situations apart internally, and the two sidecar files are what make that possible:

| situation | genotype reported | evidence needed |
|---|---|---|
| an ALT call passing the quality gate | the recorded GT | `calls` |
| an ALT call below the gate | not reported | `calls` |
| observed, and observed to be reference | `0/0` | `sites` + `regions` |
| no observation to draw on | nothing, or `./.` | `sites` + `regions` |

This is machinery, not a user-facing vocabulary: the genotype in the output already says which of these happened, so there is no state column to interpret.

### The sites file cannot be rebuilt from the calls

It is tempting to derive the site catalog by taking distinct loci out of `calls.parquet`. That works **only** when the store holds an entire joint callset. Over a subset of samples, every site where nobody in that subset carries an ALT disappears without trace — and a later query then reports "never interrogated" for a position that was in fact interrogated and found reference.

Measured on 1000 Genomes, the fraction of the true site catalog recoverable from the calls alone:

| samples in the store | catalog recovered from calls |
|---|---|
| 1,000 | 46.7% |
| 50 | 12.8% |

So at 50 samples, seven out of eight interrogated sites would silently look never-interrogated. The catalog is not redundant; it is the boundary of what the store can answer.

### Runs are not coverage

`regions.parquet` looks like a coverage track and must never be read as one. A run records that one sample was successfully called, at adequate depth, at **every catalog site in `[start, end]`**. The interval form is a *compression of a per-site fact* — "called here, and here, and here" — and it says **nothing whatsoever about the bases in between**.

This follows from what the source asserted. A plain VCF reports variants and makes no claim about any other position: an unreported base was not observed to be reference, it simply was not reported. The caller may never have looked.

Consequently a run is only meaningful at positions that appear in the site catalog, and the reconstruction checks the catalog *first and returns early* for exactly that reason. Reading a run as territory would manufacture reference observations for positions never interrogated.

Every store records `cgkit.span_semantics`:

| value | meaning |
|---|---|
| `sites` | intervals mark catalog sites only. All a plain VCF can support. |
| `blocks` | intervals came from gVCF reference blocks — positive statements about spans. Not yet produced by any converter -- `vcf-toparquet` **refuses** a gVCF rather than mis-storing one. Query a gVCF directly with `vcf-varquery`, which reads reference blocks as coverage. |

Only a `blocks` store could answer for a position absent from the catalog. Stores predating the key are read as `sites`, the conservative reading.

## What it can and cannot answer

**Can:** any question about a variant that appears in the site catalog — carriers, non-carriers, allele counts, per-sample variant lists, four-state classification.

**Cannot:** anything about a position the source VCF did not report. Such a locus yields nothing for **any** sample — not an ALT call and not a reference one — even where run intervals bracket it on both sides, and `vcf-varquery` prints a warning so that "0 carriers" is never mistaken for a real negative:

```
$ cgkit vcf-varquery --variant chr1:250:A:G --hom-ref cohort/
warning: chr1:250:A:G is not in the source; reporting not-assayed for every sample.
         A VCF only supports queries for the variants it contains.
chrom	pos	ref	alt	sample	gt	dp	min_dp	ad_ref	ad_alt	gq
```

This is a property of the *input*, not a limitation of Parquet. A gVCF, whose reference blocks carry `END` and `MIN_DP`, makes positive statements about spans -- and `vcf-varquery` now reads one directly, answering for positions no variant record mentions. What is not yet possible is *storing* that: conversion refuses a gVCF, because a block in this schema would become a `<NON_REF>` catalog site with its span discarded.

The refusal fires on **a reference block record**, not on anything in the header. A joint-genotyped cohort VCF — a DRAGEN msVCF, a GATK callset out of `GenotypeGVCFs` — routinely keeps the `##ALT=<ID=NON_REF>` declaration it inherited from its source gVCFs without containing a single block, and `##ALT=<ID=*>` declares the ordinary spanning-deletion allele. Those files convert. A record whose alternates are *all* block alleles (`<NON_REF>`, `<*>`, or a bare `.`) is what stops the conversion, and the error names the offending record. A record carrying the block allele beside a real one (`G,<NON_REF>`) is a variant record: it converts, and `<NON_REF>` is masked out of the catalog and the allele counts the way any non-focal alternate is. `-v` reports how many were masked.

### What counts as "called"

A site counts as callable for a sample only when the caller **actually made a call** there *and* depth clears `--min-dp`. Depth alone is not enough: a `./.` genotype at DP 40 means the caller saw reads and still declined, which is not the positive observation a callable region asserts. Both fixtures below exercise this:

```
chr1  400  .  T  C  .  PASS  .  GT:DP:GQ  ./.:40:0    0/0:30:99
                                          ^^^^^^^^    ^^^^^^^^^
                                          S1: DP 40   S2: called
                                          but no call
```

```
$ cgkit vcf-varquery --variant chr1:400:T:C --hom-ref --min-dp 10 cohort/
chrom	pos	ref	alt	sample	gt	dp	min_dp	ad_ref	ad_alt	gq
chr1	400	T	C	S2	0/0	.	10	.	.	.
```

### Normalization

Multiallelic records are split so each ALT allele gets its own rows. Within a split row the focal allele is recoded to `1`, reference stays `0`, and **any other alternate becomes `.`** rather than `0`. A `1/2` sample genuinely carries both alternates and must appear as a carrier in each split row; calling the other allele reference would invent a reference observation the data does not support.

Given `chr1 200 . C T,G` with `S1=1/2`, `S2=0/1`, `S3=0/0`, `S4=2/2`:

```
$ cgkit vcf-varquery --variant chr1:200:C:T cohort/
chr1	200	C	T	S1	1/.	26	26	5	12	99
chr1	200	C	T	S2	0/1	29	29	18	11	45

$ cgkit vcf-varquery --variant chr1:200:C:G cohort/
chr1	200	C	G	S1	./1	26	26	5	9	99
chr1	200	C	G	S4	1/1	27	27	0	27	99
```

`AD` is taken **per allele** (`ad_ref` is `AD[0]`, `ad_alt` is that allele's own depth), never summed: the depth supporting allele 1 says nothing about allele 2.

**Indels are not left-aligned.** Normalize beforehand if the source is not already normalized.

## AC/AN versus the sample counts

`sites.parquet` carries both, because they answer different questions and **neither can be derived from the other**.

- `ac` / `an` are **allele** counts, per the VCF convention. A `1/1` genotype contributes 2 to each. They come from `GT` alone and are not depth-gated, so `AF` is exactly `AC/AN`.
- `n_carriers` / `n_called` / `n_lowdp` are **sample** counts. A `1/1` genotype is one carrier, and `n_called`/`n_lowdp` additionally reflect the `--min-dp` used at conversion.

So `AC ≥ n_carriers` wherever a homozygote occurs, and `AN` is unrelated to `n_called` in both unit and gating. At `chr1:100 A>G` with `S1=0/1`, `S2=0/0`, `S3=1/1`, `S4=0/0`:

```
$ cgkit vcf-varquery --variant chr1:100:A:G -v cohort/
  site        AC=3 AN=8 AF=0.375  n_carriers=2 n_called=4 n_lowdp=0
  carriers    2
chrom	pos	ref	alt	sample	gt	dp	min_dp	ad_ref	ad_alt	gq
chr1	100	A	G	S1	0/1	30	30	20	10	99
chr1	100	A	G	S3	1/1	28	28	0	28	99
```

Two carriers, three alt alleles — S3 contributes two.

Both are computed **over the samples in this store**, not copied from the source's `INFO` fields, which would be wrong after splitting a multiallelic record or converting a subset of a cohort. They are also computed outside the `--no-callable` guard, since allele counts are a property of genotypes, not of coverage.

## Querying with vcf-varquery

The backend is inferred from the path — VCF or store — and the same question must give the same answer either way. That equivalence is the point of the abstraction, and it is verified on 3,202 1000 Genomes samples including the worst-covered site in BRCA1.

### Selecting sites

`--variant` is the only site selector; there is no `--region`. It takes any of these, repeatably:

```
chr1                     a whole contig
chr1:1000-2000           a region
chr1:1000                any variant at that position
chr1:1000:A:T            one exact variant
panel.vcf|.bed|.txt      a file, format detected from its content
```

A value is a file when one exists by that name and an inline selector otherwise, so a mistyped locus still gets a locus error rather than "no such file". Three file formats are recognised:

| format | shape | coordinates |
|---|---|---|
| VCF/BCF | announced by `##fileformat`; one target per ALT allele | 1-based |
| BED | `chrom start end` | **0-based half-open** |
| site list | whitespace-separated `chrom pos [ref alt]` | 1-based |

BED and site lists are told apart by their **third column**: a BED end coordinate is numeric, a REF allele never is. Both tolerate extra columns, `#` comments, blank lines and BED `track` lines, and any whitespace separates. A line holding a single token is parsed with the inline grammar above, so a file may simply list those tokens and may mix the two forms. The same file works with `vcf-gtcount --sites`, which shares the parser.

Because a misdetection would shift coordinates by one rather than fail, `-v` reports which format each file was read as and how many targets it produced, and a file yielding no targets is an error rather than an empty result.

An exact locus needs four colon-fields, a numeric position, **and** non-numeric REF/ALT. That last condition is what lets a contig name carry colons: GRCh38's ALT contigs are spelled like `HLA-A*01:01:01:01`, which also splits into four fields — but its last two are numeric. Anything not matching a locus or region shape is taken as a contig name whole, and a target that matched no rows is reported so a malformed locus does not read as a real negative. A locus *on* such a contig cannot be written inline at all; a file's columnar form takes it.

### Output layout

Every non-VCF output uses one layout, whichever selectors were named, so results can be concatenated, sorted and cut the same way:

```
chrom  pos  ref  alt  sample  gt  dp  min_dp  ad_ref  ad_alt  gq
```

The locus leads as four columns rather than one packed `chrom:pos:ref:alt` field, so rows sort and cut on position without re-splitting a composite key. By default only **valid ALT calls** are reported — genotypes carrying the alternate allele and passing the quality gate.

`min_dp` is the tightest lower bound on depth the backend can vouch for, and it exists so that a reference call still carries the evidence that the site was covered:

| row | `dp` | `min_dp` |
|---|---|---|
| a call that records its own depth | the depth | the same depth — an exact value is its own bound |
| a reference call recovered from a store | `.` | the conversion `--min-dp`, vouched for by the run it came from |
| no depth recorded anywhere | `.` | `.` — nothing is known |

The second row is the point. A store never wrote the reference genotype down, so there is no depth to report — but the sample was inside a callable run built at a known threshold, so `DP ≥ 10` is a fact the store can still assert. Writing that threshold into `dp` instead would claim a depth the data never had.

### Who carries a variant

```
$ cgkit vcf-varquery --variant chr1:500:A:T cohort/
chrom	pos	ref	alt	sample	gt	dp	min_dp	ad_ref	ad_alt	gq
chr1	500	A	T	S1	0/1	30	30	.	.	99
```

### What does a subject carry

```
$ cgkit vcf-varquery --sample S2 cohort/
chrom	pos	ref	alt	sample	gt	dp	min_dp	ad_ref	ad_alt	gq
chr1	100	A	G	S2	0/1	30	30	.	.	99
```

### Reference calls

By default a query reports only the ALT calls — "which variants does this person carry". `--hom-ref` switches it to every interrogated site — "show me all the sites for this person", reference calls included. The `gt` column separates them.

```
$ cgkit vcf-varquery --variant chr1:300:G:A --hom-ref cohort/
chrom	pos	ref	alt	sample	gt	dp	min_dp	ad_ref	ad_alt	gq
chr1	300	G	A	S1	0/0	.	10	.	.	.
chr1	300	G	A	S2	0/0	.	10	.	.	.
```

`0/0` here means the **whole genotype** was reference. At the multiallelic record above, `S2` is `0/1`: not a carrier of `G`, but not reference either, so it appears under neither. Only `S3` is genuinely reference at both split loci. Emitting `0/0` for a `0/2` sample would be a genotype the source never contained.

The dots are honest, not missing data: a store keeps only ALT genotypes, so the reference *call* was never written down — only the fact that it was made. The same query against a VCF reports the recorded genotype and its real DP/AD/GQ.

### Chromosome naming is auto-converted

UCSC (`chr22`), Ensembl (`22`) and NCBI RefSeq (`NC_000022.11`) spellings all resolve, in either direction, against a store or an index using any of them. Comparisons use canonical identity, and tabix lookups resolve against the index's own reference list. The chromosome is echoed back the way it was asked for:

```
$ cgkit vcf-varquery --variant 1:500:A:T cohort/
chrom	pos	ref	alt	sample	gt	dp	min_dp	ad_ref	ad_alt	gq
1	500	A	T	S1	0/1	30	30	.	.	99
```

Conversion preserves the source's own spelling in the store rather than rewriting it. A contig the file lacks is an absence — no rows plus a warning — but an unresolvable *region* is an error, since a region names a contig the caller asserts exists. That distinction survives `--region` being folded into `--variant`, because it keys off the selector's shape rather than a flag name: `chr1:1000:A:T` asks a question, `chr1:1000-2000` and a bare `chr1` make an assertion.

### Match the query --min-dp to the conversion --min-dp

A store baked its conversion `--min-dp` into the callable runs. Querying at a different threshold is not asking a question the store can answer consistently, and the two backends will stop agreeing on which samples are reference versus unassayed. `-v` says so rather than leaving it invisible:

```
  min-dp      10 at conversion
  NOTE: querying at --min-dp 20 but the runs were built at 10;
        non-carrier vs not-assayed will not match a direct VCF query.
```

### A gate can only act on a field the data has

`--min-dp` **fails open**: a gate over data lacking DP admits everything rather than rejecting it, because an absent depth is not evidence of a shallow one. That is deliberate, and it means a filter can silently do nothing — so both commands report field presence. At conversion:

```
fields present (a gate can only act on a field the data has)
  DP  100.0%
  GQ  100.0%
  AD  ABSENT -- per-allele depths will have no effect
```

and at query time, `--min-dp` over DP-less data warns rather than quietly returning everything.

There is deliberately **no `--min-gq`**. GQ is recorded per ALT call, but callable runs are built from depth alone, so no GQ survives for a reference call — a store could not honor the gate there while a VCF would, and the two backends would silently disagree about which samples are reference. The `gq` column is in the output, so filter on it downstream.

### Unfinished stores are unreadable

A store carries a `manifest.json.gz`, written after every member is closed. Its presence is the claim that the conversion reached the end; without it the store is refused:

```
Error: cohort/ has no readable manifest.json.gz, so it cannot be shown to be a
       complete store
       it is from an interrupted conversion, or predates manifests; re-convert
       it, or inspect it with vcf-varsummary
```

Nothing else can answer this. The parquet footers prove each member was *finished* — a footer is written only by the writer's close — but a set of finished members says nothing about how much of the input went into them. Worse, the metadata in the calls file actively misleads: `cgkit.source` and `cgkit.contigs` are stamped before the first record is read, so a store holding three chromosomes of a 22-input conversion names all 22 inputs and declares all 22 contigs. It opens cleanly, queries cleanly, and reports `not_assayed` for the rest — exactly what a complete store reports for a position the source never mentioned.

So the manifest also records what was *written*: per-member row counts, checked at open against each member's own footer, and a per-chromosome census of sites and calls. That census is the one field that can contradict the rest, and `vcf-varsummary --counts` prints it.

There is no escape hatch, and a store written before manifests must be re-converted.

Note this is a different thing from the section below, which is about a store that was finished but tracked no coverage.

### Stores that cannot classify refuse rather than degrade

If `sites.parquet` or `regions.parquet` is missing, or the store was built with `--no-callable`, then `--hom-ref` returns `ErrNotClassifiable` instead of an answer:

```
Error: store cannot distinguish non-carrier from not-assayed: cohort/regions.parquet is missing
```

Reporting an unobserved sample as reference would invent an observation, so the failure is loud.

## Performance

### Row-group pruning

Input is coordinate-sorted and written in stream order, and the writer declares that ordering, so parquet's per-row-group `min`/`max` on `pos` are tight and nearly non-overlapping. A locus lookup reads the footer, skips every row group whose position range cannot contain the target, and decodes only the survivors.

The predicates are conservative in exactly one direction: **they may keep a row group that turns out to hold nothing, but they must never skip one that holds a match.** The per-row checks downstream are unchanged, so pruning can only affect how much is read, never what is found.

Chromosome bounds are used **only** when a row group holds exactly one chromosome, because the statistics are lexicographic while chromosome names are not ordered that way — comparing `chr17` against a `chr1`..`chr9` range would be meaningless. Restricting the test to single-chromosome groups keeps it sound and still catches the common case, since coordinate-sorted input puts almost every group inside one chromosome.

**Supply multi-VCF input in coordinate order.** Correctness does not depend on it — the answers are identical either way — but out-of-order input widens the per-group position ranges and pruning gets weaker.

### Historical measurements

Every number in this table is an ad-hoc measurement recorded in a commit message, taken on 1000 Genomes data that is **not committed to this repository**, so none of it can be re-checked. The sections after it come from `go test -bench` over a generated corpus and can be. Read the middle column carefully; these measure different things and are easy to conflate.

| result | what was measured | source |
|---|---|---|
| 479ms → 157ms (**3.05x**) | one locus lookup, before vs after row-group pruning, on a two-chromosome store of 819,389 calls | commit `7f4ca94` |
| 298ms → 166ms (**~1.8x**) | one locus lookup, inputs supplied out of coordinate order vs in order | commit `1aa5ce2` |
| 429ms vs 426ms (**no effect**) | a `--sample` query with vs without the `sample_id` bloom filter | commit `7f4ca94` |
| 417ms vs pyarrow's 26.5ms | reading all columns of the same store — i.e. Go row-by-row decode overhead | ad-hoc, same store |

### Store versus plain VCF

Measured, and the answer is more qualified than the direction alone suggests. Reproduce with:

```
CGKIT_BENCH_SAMPLES=100 CGKIT_BENCH_SITES=5000 \
  go test ./internal/cmd/vcfcmd/ -run '^$' -bench BenchmarkStoreVsVcf
```

100 samples × 5,000 sites, ~5% of genotypes carrying the alternate, so ~25,000 calls in the store:

| targets | store | VCF, unindexed | VCF, tabix-indexed |
|---|---|---|---|
| 1 locus | **7.1 ms** | 310 ms | **13.3 ms** |
| 100 loci | **17.1 ms** | 325 ms | 312 ms |

Two things worth reading carefully.

**Against an indexed VCF a single-locus lookup is only about 1.9× faster, not 40×.** The 310 ms unindexed figure measures a file that has to be read end to end because nothing tells the reader where to look — that is a property of the file not being indexed, not of the format. Quoting it as the comparison would flatter the store considerably.

**The 18× at 100 loci is partly a limitation of this implementation, not of the formats.** `VcfStore` only seeks when a query names exactly one locus or one region; for several it scans and filters, so the index goes unused. A per-locus seeking path would narrow that gap. The store's advantage that *is* structural is holding ~4% as many genotypes and reading only the columns named.

### Bulk queries versus one query per locus

This is the measurement that justifies `Calls` taking a whole target set rather than being called in a loop:

```
go test ./internal/cmd/vcfcmd/ -run '^$' -bench BenchmarkQueryPanel
```

| targets | one bulk query | one query per locus | ratio |
|---|---|---|---|
| 1 | 7.0 ms | 7.3 ms | 1× |
| 10 | 15.1 ms | 67.7 ms | 4.5× |
| 100 | 17.0 ms | 662 ms | 39× |
| 1000 | 18.9 ms | 6,612 ms | **349×** |

The bulk path is nearly flat — 7 ms to 19 ms across a thousandfold increase in targets — because it makes one ordered pass whatever the panel size. The loop is linear in targets, since each lookup re-opens the store and re-parses the footer. The crossover is somewhere between 1 and 10 targets, so looping loses almost immediately.

Extrapolating the loop to a 10⁶-variant panel gives tens of hours, which is why the API does not expose a shape that invites it.

### Reading these numbers

The corpus is **synthetic and small**, and the machine is a QEMU virtual CPU, so the absolute milliseconds mean little — the ratios and the *shapes* (flat versus linear) are the findings. Row-group size matters enough to be a knob: left at the converter's 250,000 default, a corpus this size lands in a single row group, position statistics can never exclude it, and every lookup scans the whole file. `CGKIT_BENCH_ROWGROUP` controls it, and the benchmark defaults to a value that produces several groups, as a real store would have.

Deliberately not reported: bytes read and row groups decoded. Those would be the more durable numbers, since they survive hardware changes, but they need counters inside `varstore`.

### Known limitation: --sample queries do not prune

Every sample appears in nearly every row group, so the `sample_id` bloom filter always answers "maybe" and no group can be skipped. Making these fast would need a sample-major sort, which directly conflicts with the `(chrom, pos)` sort that locus pruning depends on.

Measurement says the layout is not even the main cost: the same store reads in 26.5ms under pyarrow versus 417ms in Go, so roughly 16x of it is row-by-row struct materialization, not storage. A sample-sorted second copy would double the largest file to fix a problem that mostly is not storage, so it is deliberately deferred in favour of two zero-storage fixes — column projection (read `sample_id` alone to locate rows, 4.3x less data) and the decode path itself.

## Safety properties

Two failure modes are guarded structurally rather than by care.

**Conversion refuses to overwrite an existing store.** If any member is already present under `--out`, conversion stops and asks for `--force`. Writing truncates them all, and a half-replaced set is worse than either keeping or replacing the old one. The guard keys on the members, so an existing directory holding unrelated files is fine and its contents are left alone.

**A conversion that cannot finish must not leave something that looks finished.** Every error path calls `Writer.Discard()`, which unlinks the members *without* finalizing them — closing a parquet writer writes a complete, valid footer, so finalizing files that are about to be removed meant a process killed mid-discard left behind exactly the well-formed partial store this prevents. `Close` likewise stops at the first failure instead of finishing the other members: three structurally valid files of which one is silently short is the one outcome a reader cannot detect. A failed removal is reported, because what survives is what blocks the retry.

## Multiple inputs

Whole-genome callsets usually ship as one VCF per chromosome, so several inputs may be given and are concatenated into one store.

They must carry **exactly the same samples**. Differing column order is **remapped**, not merely tolerated: genotype columns are addressed positionally, so getting this wrong does not error — it silently attributes every genotype to the wrong person and produces entirely plausible output. A sample-set mismatch is an error naming what differs.

Inputs must **not overlap**: a chromosome cannot be revisited once left, and positions cannot go backwards within one. Overlapping inputs would write the same site twice and split its AC/AN across two rows.

## Flags

### vcf-toparquet

| flag | default | description |
|---|---|---|
| `--out DIR` | *required* | the store directory, created if needed |
| `--force` | off | overwrite an existing store at `--out` |
| `--min-dp N` | `10` | depth at or above which a site counts as callable |
| `--no-callable` | off | accept a source with no DP field; regions will be empty |
| `--passing` | off | skip filtered records |
| `--region R` | | only variants in this 1-based region; needs a tabix index |
| `--compression C` | `zstd` | `zstd`, `snappy`, or `none` |
| `--row-group-size N` | `250000` | rows per parquet row group |
| `-v`, `--verbose` | off | progress, per-chromosome trace, conversion summary, field census |

### vcf-varquery

| flag | default | description |
|---|---|---|
| `--variant TARGET` | | a locus, region, contig, or a file of them (repeatable) |
| `--sample NAME` | | report variants carried by this subject (repeatable) |
| `--min-dp N` | `0` | minimum DP for a call to count |
| `--hom-ref` | off | report every interrogated site, not only the ALT calls |
| `--tbi` | off | also write a tabix index (needs `-o` with a `.gz` name) |
| `--format F` | `tsv` | `tsv`, `json`, `vcf` (sample mode), or `list` (variant mode) |
| `--store KIND` | infer | force the backend: `vcf` or `parquet` |
| `-o`, `--output` | `-` | output filename |
| `-v`, `--verbose` | off | backend, provenance, per-variant AC/AN/AF, gate effect (all to stderr) |

### vcf-varsummary

| flag | default | description |
|---|---|---|
| `--samples` | off | list the sample roster, one per line |
| `--sites` | off | stream the variant catalog (a full pass for a VCF) |
| `--counts` | off | totals and the per-chromosome census (a full pass for a VCF) |
| `--format F` | `text` | `text`, or `json` to emit a store's manifest verbatim |
| `--store KIND` | infer | force the backend: `vcf` or `parquet` |
| `-o`, `--output` | `-` | output filename |
| `-v`, `--verbose` | off | note what is being read, and what a scan will cost |

The default report reads no records at all — it is O(header) on both backends, so it is instant against a whole-genome VCF. Everything requiring a pass over the data is opt-in. A store answers from its manifest; a VCF has none, since it is not a conversion, and the report says so rather than inventing a `--min-dp` it cannot know. For an indexed VCF the contigs that carry records come from the tabix index for free.

All verbose output goes to stderr, so the tabular stream stays parseable.

## Provenance

Conversion records these keys in the calls file, reported by `vcf-varquery -v`:

| key | contents |
|---|---|
| `cgkit.samples` | the sample roster, in source order |
| `cgkit.min_dp` | the callable threshold used |
| `cgkit.span_semantics` | `sites` or `blocks` |
| `cgkit.contigs` | the source's `##contig` lines, verbatim — so an export can declare its reference |
| `cgkit.nocallable` | `1` when regions are absent by request |
| `cgkit.source` | input filename(s) |
| `cgkit.program`, `cgkit.command` | cgkit version and full command line |

The roster lives in metadata rather than a fourth file because reconstructing reference calls needs every sample, while the calls file only ever mentions carriers — a sample carrying nothing anywhere would otherwise be invisible and never reported at all.

## See also

- `cgkit vcf-toparquet --help`, `cgkit vcf-varquery --help`
- `cgkit vcf-gtcount` — genotype distribution summarized per site
- The implementation lives in the `varstore` package; its doc comments carry the design rationale in more detail.
