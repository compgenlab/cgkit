# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

cgkit is a Go CLI toolkit for computational genomics research. It provides commands for sequence analysis, NGS data wrangling, and bioinformatics operations. The underlying library (sequence I/O, alignment, SAM/BAM/CRAM handling) lives in the separate `cghts` module (`github.com/compgenlab/cghts`).

**The genotype store is not here.** `varstore` (`github.com/compgenlab/varstore`) is its own binary for creating, manipulating and querying varstores: `vcf-tovarstore`, `vcf-varstore create`, `vcf-refblocks`, `vcf-varstore attach-refblocks`, `vcf-varquery` and `vcf-varsummary` moved there (as `convert`, `create`, `refblock-calc`, `refblock-attach`, `query`, `summary`) and were removed from cgkit; do not re-add them. The store format itself is `cghts/varstore`. `vcf-gtcount` stays here and shares the store's site-selector grammar through `varstore.ParseTargets`.

**Oxford Nanopore UMI and poly(A) tooling is not here.** It lives in `nupa` (`github.com/compgenlab/nupa`, private), built on the same library. The `ont-*` commands were removed from cgkit; do not re-add them. `ont-umi-lookup` was withdrawn rather than migrated and exists in neither repo.

**Module:** `github.com/compgenlab/cgkit`
**Go version:** 1.24.9
**CLI framework:** Cobra
**Library dependency:** `github.com/compgenlab/cghts`

## Commands

```bash
# Build all targets (darwin_arm64, darwin_amd64, linux_arm64, linux_amd64, windows_amd64)
make build

# Run all tests
make test
# equivalent to:
GOCACHE=/tmp/go-build-cache go test ./...

# Run a single test
go test ./internal/cmd/samcmd/... -run TestSamStats

# Run with CPU profiling
./cgkit --profile=cpu.prof <subcommand>
```

## Dependency on cghts

All format I/O and algorithms come from `github.com/compgenlab/cghts` (packages
`seqio`, `align`, `htsio` and its subpackages, `support/*`, `analysis/seq`).

How the `cghts` dependency resolves, by context:
- **Local builds** use the `go.work` workspace (parent dir, untracked) that joins
  a sibling `cghts` checkout, so you build against your live local `cghts` tree. The
  `Makefile` deliberately does **not** set `GOWORK=off`.
- **Remote/CI builds** (no `go.work` present) use the **latest released cghts from
  GitHub**: the GitHub Actions workflow runs `go get github.com/compgenlab/cghts@latest`
  before vet/test/build, with `GOPRIVATE=github.com/compgenlab/*` so a freshly
  pushed cghts tag is fetched directly from GitHub (no module-proxy/sumdb lag).
- The committed `go.mod` pin is the fallback for `go install` users and the
  source archive; keep it current with `make bump-cghts`.
- **Building against unreleased cghts.** A feature spanning both repos can pin a
  cghts *commit* instead of a tag:
  `GOWORK=off go get github.com/compgenlab/cghts@<sha> && GOWORK=off go mod tidy`,
  which writes a pseudo-version (`v0.12.1-0.<date>-<sha>`). CI honours any pin
  carrying a `-` (pseudo-version or pre-release) rather than following `@latest`,
  which skips those and would build the branch against code lacking what it needs;
  a release pin (`vX.Y.Z`) still means `@latest`. Re-pin as the cghts side moves.
  `make release-check` **refuses** a pseudo-version pin, so tag cghts and
  `make bump-cghts` before releasing.

### Cutting a release
The cghts tag must land on GitHub before cgkit builds against it:
1. **cghts**: tag `vX.Y.Z` on `main`, push the tag.
2. **cgkit**: `make bump-cghts` (pins `go.mod` to the new cghts), commit
   `go.mod`/`go.sum`, push. CI's `go get cghts@latest` then resolves the same tag.

## Architecture

This repo holds only the CLI layer: `main.go` (entry point with `--profile`
support) and `internal/cmd/` (Cobra commands). The only third-party dependencies
are cobra/pflag; everything genomics-related is delegated to `cghts`.

### Remote inputs

Inputs may be local paths, `http(s)://` URLs or `s3://` objects; **outputs are
local** (writing a varstore to S3 is `varstore convert`'s, not cgkit's). Dispatch is cghts's `iosource`, and index sidecars resolve over
the same transport as the data, so an indexed region query still seeks.

- The `s3://` blank import lives in **`internal/locator`**, a leaf package, not in
  `main.go`. `main` is not linked into `go test ./internal/...`, so a test
  asserting the "no transport registered (registered: s3)" message would pass
  there while the shipped binary said otherwise; keeping it in a package the
  tests link makes the test the guard on the import. This is why cgkit carries
  the AWS SDK: measured **26.5 MB → 38.4 MB** for linux_amd64 (+11.9 MB, +45%).
  `-ldflags "-s -w"` brings that back to 26.8 MB — not enabled, since stripping
  is an orthogonal decision, but the option is there if size matters.
- `internal/locator.IsRemote` wraps `iosource.IsRemote` rather than
  reimplementing it. cgkit must classify a locator exactly as the library will
  dispatch it; a classifier that disagreed by one edge case (a Windows drive
  letter is not a scheme) would not fail loudly, it would open the wrong thing.
- **Site selectors classify scheme-first, then `os.Stat`** (`varstore.IsTargetFile` in cghts, used by `vcf-gtcount`). A URL cannot be
  stat'd, and a stat-first rule would demote it to an inline selector, which
  parses as a contig name and matches nothing — an empty result rather than an
  error. Only values containing `://` past the first character change meaning,
  and nothing in the selector grammar can look like that, so a mistyped locus
  still gets a locus error and `HLA-A*01:01:01:01` is untouched. Both pinned by
  tests; this is the single riskiest change in the remote work.
- **Outputs are refused by name**, via `locator.CheckLocalOutput`, always *after*
  a command's stdout branch or `-` breaks. A BAM or VCF writer has no
  writer-side equivalent of a byte source.
- Costs stated rather than hidden: a *streaming* read of a remote object
  transfers the whole thing (no index to skip with), and `vcf-gtcount --threads N`
  fetches the index N times because a tabix reader is not concurrency-safe.
- Deliberately **not** covered: FASTA/FASTQ/BED streaming inputs, `vcf-concat`'s
  chunked reader, `vcf-annotate` sources, remote writes, `gs://`.

### CLI Command Structure

Commands are registered in `internal/cmd/root.go` and grouped by file format or domain:
- `fasta-cat`, `fasta-wrap`, `fasta-gc` — FASTA operations
- `fastq-gc` — FASTQ operations
- `sam-stats` — Summary statistics for SAM/BAM/CRAM: read counts, mapping rates, Q30, depth, SAM flag breakdown, per-reference counts, optional `--tags` value distributions and `--calc-insert` median. Profiles the first read of each pair only (ports `ngsutils bam-stats`). Phase 1 omits the `--gtf` gene-model and `--bed` on-target stats.
- `sam-vcfcount <bam> <vcf>` — Per-variant ref/alt read counts from an indexed SAM/BAM/CRAM, one row per ALT (ports ngsutilsj `vcf-count`, but with the alignment file first, like every `sam-*` command). **Reads are classified whole, not per pileup column**: the read's sequence over the REF span (aligned bases plus any insertion inside or right after it) is compared to REF and each ALT, so indels, MNVs and complex alleles match exactly instead of through mpileup's `+N`/`-N` encoding. An indel's read must reach the base after REF (else a read ending on an insertion's anchor reads as REF), and a read whose deletion covers a flank is *other*, not ALT (else a longer deletion matches a shorter one). A fragment counts once when both mates overlap; mates that disagree count for neither (`--both-mates` disables). `--min-basequal` applies to SNVs/MNVs only, as in mpileup. `--flag-exclude` defaults to 3844 (unmapped/secondary/QC-fail/dup/supplementary). No realignment: an indel the aligner placed elsewhere in a repeat is missed. Contig names resolve through `ContigConverter`; a contig absent from the BAM is an error unless `--skip-missing`. Symbolic alleles, breakends, `*` and gVCF blocks are not counted. Reads are fetched once per `--batch-len` window, and the read loop relies on query results being coordinate-sorted. `--min-depth N` writes only variants with `total_count` ≥ N (1 drops the uncovered `0 0` rows); depth is every spanning read, *other* included, so a site covered only by a third allele counts as covered, and it is decided per record so a multiallelic record's rows are kept or dropped together. `--vcf-counts`/`--vcf-af` add the caller's own view from the `--sample` `AD` (ngsutilsj column order: after `ref alt`, before the read counts); `vcf_alt_freq` is alt/(ref+alt) from AD, unlike `alt_freq`'s spanning-read denominator. FORMAT/AD must be declared (refused up front); a record lacking AD reports `.` rather than failing as ngsutilsj did. `--snvs` is strict (single-base REF *and* ALT, so MNVs are out) and, like `--indels`, filters **per allele, not per record**: `A>G,AT` keeps the G row. The filter only hides rows -- a hidden ALT is still matched, because unmatching it would let a read carrying that insertion fall through to REF and inflate `ref_count` on the rows that remain.
- `seq-pairwise`, `seq-revcomp` — Sequence analysis
- `vcf-split` — Split a VCF into `BASE.N.vcf.gz` chunks of N variants each. **A series is written whole or not at all**: it refuses an existing series without `--force`, `--force` clears the *entire* old series rather than only the chunks the new run overwrites, and a run that fails partway removes what it wrote. All three exist because `vcf-concat --chunks` stops at the first *missing* file rather than the first invalid one — so a truncated tail chunk, or a leftover chunk from a longer previous run, is silently reassembled as though it were the whole input. A refused `--tbi` index does **not** abort the split: `writeVcfTbi` rejects an unsorted chunk *after* that chunk is complete, and stopping there would leave exactly the partial series the guards prevent, so the series is finished, the indexes are skipped, and the error is reported at the end.
- `vcf-strip` — Remove annotations (FILTER/INFO/FORMAT/samples/ID) from a VCF. **Refuses to drop `INFO/END` from a gVCF** unless `--force-end`, because a reference block's extent lives entirely in that field: removing it leaves a file that still parses while claiming one base of coverage where it claimed thousands. Detection is from the header (`##ALT=<ID=NON_REF>` or `<*>`, a `##GVCFBlock` line, a declared `MIN_DP`/`LEN`) because the output header is written before the first record is read, so a record-level answer arrives too late to keep the `##INFO` declaration; if detection misses and a block carrying `END` turns up anyway, it fails rather than emit records with an undeclared `END`. `--all --keep-info END` was always an escape hatch; `--force-end` is the inverse.
- `vcf-tstv` / `vcf-stats` — Both report `-` for a Ts/Tv ratio with no transversions to divide by. This **diverges from ngsutilsj deliberately**, which is not self-consistent here: its `vcf-stats` emits an empty field and its `vcf-tstv` emits `Infinity` for the same input, one reading as a failure and the other as a real number. `testdata/sample.vcf` has no transversions, so every parity case hits this; the parity harness normalizes the reference output for that one field rather than pretending the divergence is absent.
- `vcf-tobed` — Export variant positions as BED. The interval is the **reference extent** the record covers (`rec.RefSpan()`), not the alternate position: `AltPositions` resolves an SV's *partner breakpoint*, which may be on another chromosome, and using it as an extent gave a plain deletion one base too many and collapsed a multi-base-ALT deletion to a single base. Reference blocks are skipped -- they describe coverage, not variants -- and a mixed `G,<NON_REF>` record keeps its real allele. This **diverges from ngsutilsj deliberately**, which ends a deletion at `POS+len(REF)`; the parity goldens still pass only because `testdata/sample.vcf` has no multi-base-REF deletion.
- `vcf-gtcount` — Genotype (GT) distribution per site across samples, for sites given inline or via `--sites`. Answers a different question from `varstore query` (a per-site *summary*, not per-sample genotypes) and looks each site up by random access, so it needs a tabix-indexed input and **refuses** a region, whole contig or BED interval — there is no single position to look up. Shares `varstore query`'s target grammar and parser (`varstore.ParseTargets` in cghts).
