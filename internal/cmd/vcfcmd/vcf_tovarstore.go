package vcfcmd

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/compgenlab/cghts/varstore"
	// Registers the s3 sink, so --out may name a bucket. A blank import:
	// nothing here calls it, and without it s3:// is refused as an unknown
	// scheme rather than written.
	_ "github.com/compgenlab/cghts/varstore/sinks3"
	"github.com/compgenlab/cghts/vcf"
	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/spf13/cobra"
)

var (
	vcfToVarstoreOut          string
	vcfToVarstoreMinDP        int
	vcfToVarstoreNoCallable   bool
	vcfToVarstoreCompression  string
	vcfToVarstoreRowGroupSize int
	vcfToVarstoreForce        bool
	vcfToVarstoreInfo         []string
	vcfToVarstoreBands        []int
	vcfToVarstoreShardSites   int
	vcfToVarstoreFormat       []string

	// One string per reserved key, keyed by the key itself, plus the generic
	// --meta pairs. Both feed one map; see collectMeta.
	vcfToVarstoreMetaNamed = map[string]*string{}
	vcfToVarstoreMetaPairs []string
)

// collectMeta folds the named --meta-<key> flags and the generic --meta
// KEY=VALUE pairs into the single map the store records.
//
// A key given twice is an error rather than last-wins, including across the two
// spellings (--meta-dataset X --meta dataset=Y). A store's metadata is a claim
// about what it holds, and silently picking one of two conflicting claims is
// the class of thing this command already refuses elsewhere -- it errors on a
// sample-set mismatch rather than guessing which roster was meant.
func collectMeta() (map[string]string, error) {
	meta := map[string]string{}
	from := map[string]string{} // key -> the flag that set it, for the error

	set := func(key, val, via string) error {
		// Checked here as well as by the writer, which validates for every
		// caller. Doing it in the pre-flight means the run dies before
		// EnsureStoreDir, so a typo leaves nothing behind, and the message can
		// name the flag that carried the key -- which the library cannot.
		if !varstore.ValidMetaKey(key) {
			return fmt.Errorf("%s: invalid metadata key %q: keys must be non-empty and use "+
				"only lowercase letters, digits, underscore and hyphen", via, key)
		}
		if prev, ok := from[key]; ok {
			return fmt.Errorf("metadata key %q given twice, by %s and %s: "+
				"remove one rather than leaving which value is recorded to chance", key, prev, via)
		}
		from[key] = via
		meta[key] = val
		return nil
	}

	for _, key := range varstore.ReservedMetaKeys {
		if p := vcfToVarstoreMetaNamed[key]; p != nil && *p != "" {
			if err := set(key, *p, "--meta-"+key); err != nil {
				return nil, err
			}
		}
	}
	for _, pair := range vcfToVarstoreMetaPairs {
		key, val, ok := strings.Cut(pair, "=")
		if !ok {
			return nil, fmt.Errorf("--meta %q: expected KEY=VALUE", pair)
		}
		if key == "" {
			return nil, fmt.Errorf("--meta %q: empty key", pair)
		}
		if err := set(key, val, fmt.Sprintf("--meta %s=...", key)); err != nil {
			return nil, err
		}
	}
	if len(meta) == 0 {
		// Nil, not an empty map: the store omits absent metadata entirely, so
		// that "not stated" stays distinguishable from "stated as nothing".
		return nil, nil
	}
	return meta, nil
}

var vcfToVarstoreCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": "v0.5.0"},
	Use:         "vcf-tovarstore <input.vcf> [input2.vcf ...]",
	Short:       "Convert a VCF to a varstore",
	Long: `Convert a VCF into a VARSTORE: a columnar genotype store that keeps only the
alternate-allele calls, along with enough context to still tell a
confidently-called reference apart from a position that was never assayed.

Several inputs may be given, which is how whole-genome callsets usually ship --
one VCF per chromosome. They must carry exactly the same samples; differing
column order is remapped, since genotype columns are positional and getting that
wrong would silently attribute every genotype to the wrong person. A sample-set
mismatch is an error naming what differs.

Inputs must not overlap: a chromosome cannot be revisited once left, and
positions cannot go backwards within one. Overlapping inputs would write the
same site twice and split its AC/AN across two rows, so this is refused.

Give them in coordinate order for best query performance. Correctness does not
depend on it -- the answers are identical either way -- but the per-row-group
position statistics stay tighter, and a locus lookup then skips more of the
file. Measured on a two-chromosome store, supplying them out of order cost
about 1.8x on a locus query (166ms against 298ms).

The name is the STORE rather than its encoding. Parquet is how a varstore is
written today, behind an interface this package owns, and the design keeps room
for another backend -- so naming the command after the format would be a promise
about something deliberately left changeable. Several varstores, disjoint by
chromosome, can then be grouped into a VARSET and queried as one -- see
cgkit vcf-varset.

The store is written to --out, and its members form one inseparable set:

  cohort/
    calls.parquet      one row per ALT-carrying genotype
    sites.parquet      one row per interrogated site, with AC/AN and counts
    regions.parquet    contiguous runs of adequately-covered sites, per sample
    manifest.json.gz   written last; the store is unreadable without it

A store is a directory, created if needed; a trailing "/" on --out is optional
and means nothing. The members are only meaningful together, so this keeps the
set as one thing to copy, move or delete, and it is the layout every Parquet
tool expects -- DuckDB or pyarrow can be pointed straight at a member.
vcf-varquery accepts the directory, with or without the slash, or any member
path within it.

The manifest is what makes a store readable rather than merely present. It is
written after every member is closed, so its presence means the conversion
reached the end -- which nothing else can tell you. The parquet footers prove
each member was finished, but a set of finished members says nothing about how
much of the input went into them, and a store that covered three of twenty-two
chromosomes answers "not assayed" for the rest, exactly as a complete store
answers for a position the source never reported. So the manifest also records
what was written: per-chromosome site and call counts, per-member row counts,
the sample roster. Read it with vcf-varsummary.

A store written by an older cgkit has no manifest and must be re-converted.

Conversion refuses to overwrite an existing store: if any member is already
present under --out it stops and asks for --force. Writing truncates them all,
and a half-replaced set is worse than either keeping or replacing the old one.
The check keys on the members, so an existing directory holding unrelated files
is a fine target and its contents are left alone.

The sites file carries both allele counts (AC, AN) and sample counts
(n_carriers, n_called, n_lowdp). They are not interchangeable: a 1/1 genotype is
one carrier but two alt alleles, so AC >= n_carriers wherever a homozygote
occurs, and AN counts alleles without regard to depth while n_called counts
samples that cleared --min-dp. AF is exactly AC/AN. Both are computed over the
samples in this store rather than copied from the source's INFO fields, which
would be wrong after splitting a multiallelic record or converting a subset of a
cohort.

The sites file is not redundant with the calls. Deriving the site list from the
distinct loci in the calls only works when the store holds an entire joint
callset; over a subset of samples, every site where nobody in that subset
carries an ALT disappears, and a later query would report those positions as
never interrogated rather than as observed and reference.

Records are normalized to one variant per row: a multiallelic record is split
so each ALT allele gets its own rows. Within a split row the focal allele is
recoded to 1, reference stays 0, and any other alternate allele becomes "." --
so a 1/2 sample is correctly a carrier of both alleles. AD is taken per allele
(ad_ref is AD[0], ad_alt is that allele's own depth) rather than summed, since
depth supporting one alternate says nothing about another. Indels are NOT
left-aligned; normalize beforehand if the source is not already.

The regions file records, per sample, runs of catalog sites at which that
sample was successfully called at DP >= --min-dp. The interval form is only a
compression of that per-site fact; it makes NO claim about the bases between
those sites.

This bounds what the store can answer. A plain VCF reports variants and says
nothing whatsoever about any other position -- an unreported base was not
observed to be reference, it was simply never reported. The sites catalog is
therefore the exact boundary of what is knowable, and a query for a locus
outside it returns not-assayed for every sample rather than a set of reference
calls, even where run intervals appear to bracket it. Only a gVCF, whose
reference blocks carry END and MIN_DP, makes positive statements about spans and
could answer off-catalog positions.

  --out DIR             the store directory, created if needed (required).
                        May be s3://bucket/store, which streams the members
                        straight there and needs no local scratch.
  --force               overwrite an existing store at --out
  --min-dp N            depth at or above which a site counts as callable
  --info ID[,ID...]     capture these INFO fields into the sites catalog
  --no-callable         proceed when the input has no DP field at all
  --passing             skip filtered records
  --compression C       zstd (default), snappy, or none
  --row-group-size N    rows per parquet row group
  -v, --verbose         report progress and a conversion summary on stderr,
                        including which of DP/GQ/AD the input actually carries

--info captures source INFO fields into sites.parquet, each as its own typed
column named info_<id> in lower case. It is for the things that are properties
of THIS FILE rather than of the variant -- imputation quality (R2, DR2, INFO),
an IMP or TYPED flag, a panel allele frequency, VQSLOD. Those cannot be
recovered from an annotation source later, because a different VCF of the same
cohort gives different numbers at the same locus, and without capture they are
simply dropped at conversion.

  --info R2,IMP         two fields by name
  --info 'AF_*'         a glob over the header's declared fields
  --info '*'            everything the header declares that can be stored

Types and cardinality are read from the ##INFO header lines, never from the
command line: the file already declares them, so restating them would only be a
way to get them wrong. Only Number=1 (one value per site), Number=A (one per
ALT) and Flag can be captured -- a site row is one ALT, so Number=R and Number=G
have values with nowhere to go. Naming such a field is an error; one that merely
matched a glob is skipped with a note.

Values are stored per ALT where the field is Number=A, so a multiallelic record
that splits into three site rows gives each row its own value rather than
repeating the first. A field absent from a record is stored as null, which stays
distinct from a value of zero -- "no R2 was reported here" and "R2 is 0 here"
are different claims.

Captured columns never collide with the store's own: --info AC lands in info_ac
and leaves the computed ac alone. That matters because AC/AN/n_called are
recomputed over the samples in the store rather than copied from the source,
which is the only reason they survive a subset.

The manifest records which fields were captured, with their source key, type and
Number. That is the only place absence is answerable: a reader gets zero both
for a column that is not in the file and for one holding zero.

Metadata records what the store *is*, as opposed to how it was made. The store
already knows the latter -- the command, the inputs, the sample roster, when it
ran -- but nothing about the former. A whole-genome callset arrives as one VCF
per chromosome and converts into one store, which can then name its 24 input
filenames but not the release they collectively are; and those filenames stop
identifying anything once the store moves. The assembly is the same gap: the
##contig lines carry lengths that imply GRCh37 or GRCh38, but a store read
against the wrong assumption does not fail, it answers with coordinates that
mean something else.

  --meta-dataset NAME       the release or callset this store was built from
  --meta-reference NAME     assembly the calls were made against
  --meta-caller NAME        variant caller and version
  --meta-accession ID       study or dataset accession
  --meta-url URL            where the source data came from
  --meta-version V          the dataset's own release version
  --meta-description TEXT   free text
  --meta KEY=VALUE          anything else (repeatable)

Values are recorded verbatim -- cgkit cannot know whether "GRCh38" is true, and
normalizing it would turn your claim into its own. Keys are lowercase
[a-z0-9_-]. A key given twice, by either spelling, is an error rather than
last-wins: which of two conflicting claims gets recorded should not depend on
the order the flags were typed. vcf-varsummary prints all of it, and
--format json emits it verbatim for jq.`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		if vcfToVarstoreOut == "" {
			return fmt.Errorf("you must specify an output store directory with --out")
		}
		// A store may be written to an object store as well as to a directory,
		// so --out is checked against the transports linked in rather than
		// refused for being remote. An unregistered scheme is named as such:
		// left to the writer it would surface as a not-found against a path
		// nobody typed.
		if !varstore.CanWrite(vcfToVarstoreOut) {
			return fmt.Errorf("--out: cannot write a store to %q; this build writes to a directory%s",
				vcfToVarstoreOut, writableSchemes())
		}
		if vcfToVarstoreMinDP < 0 {
			return fmt.Errorf("--min-dp must not be negative")
		}
		if vcfToVarstoreRowGroupSize <= 0 {
			return fmt.Errorf("--row-group-size must be a positive number")
		}
		codec, err := varstore.CodecFor(vcfToVarstoreCompression)
		if err != nil {
			return err
		}
		// Resolved before any input is opened: a conflicting or malformed --meta
		// is a typo on the command line, and there is no reason to read a
		// whole-genome callset before reporting one.
		meta, err := collectMeta()
		if err != nil {
			return err
		}

		// The first input fixes the sample roster; every later one must carry
		// the same people, though not necessarily in the same column order.
		first, err := openRecordSource(cmd, args[0], vcfRegion)
		if err != nil {
			return err
		}
		samples := first.header.Samples()
		// Resolved from the header while it is still open, and BEFORE the
		// overwrite check below: a misspelled --info must not be the thing that
		// discovers itself after the previous store has been truncated.
		infoFields, infoSkipped, err := resolveInfoFields(first.header, vcfToVarstoreInfo)
		if err != nil {
			first.close()
			return err
		}
		formatFields, formatSkipped, err := resolveFormatFields(first.header, vcfToVarstoreFormat)
		if err != nil {
			first.close()
			return err
		}
		first.close()
		// Never silently: a glob that quietly dropped a field would leave the
		// caller looking for a column that was never written.
		if len(formatSkipped) > 0 {
			fmt.Fprintf(cmd.ErrOrStderr(),
				"note: not capturing %s -- only Number=1 and Number=A fit a call row\n",
				strings.Join(formatSkipped, ", "))
		}
		if len(infoSkipped) > 0 {
			fmt.Fprintf(cmd.ErrOrStderr(),
				"note: not capturing %s -- only Number=1 and Number=A fit a site row\n",
				strings.Join(infoSkipped, ", "))
		}
		if len(samples) == 0 {
			return fmt.Errorf("%s has no samples; a genotype store needs per-sample calls", args[0])
		}
		// A gVCF is refused, but on the evidence of a reference block record rather
		// than of the header -- see gvcfRefBlockError in gvcf.go.

		// The destination, opened once and reused: the overwrite check and the
		// writer must agree about where the store is going, and opening it
		// twice is how they would come to differ.
		sink, err := varstore.OpenSink(vcfToVarstoreOut)
		if err != nil {
			return err
		}
		// Refuse to clobber an existing store before opening anything: the
		// writer truncates every member, so this is the last moment the
		// previous one still exists.
		//
		// Asked of the SINK the writer will use, not of the base again: opening
		// the destination twice is how the check and the writer come to
		// disagree about where the store is going.
		if err := varstore.CheckStoreTargetIn(sink, vcfToVarstoreForce); err != nil {
			return err
		}

		// Bands must be ascending and above the gate: a boundary below --min-dp
		// names a class nothing can fall into, which is not an error so much as
		// a sign the caller meant something else.
		bands := make([]int32, 0, len(vcfToVarstoreBands))
		for i, b := range vcfToVarstoreBands {
			if b <= 0 {
				return fmt.Errorf("--depth-bands: %d is not a depth", b)
			}
			if i > 0 && int32(b) <= bands[i-1] {
				return fmt.Errorf("--depth-bands must ascend; %d follows %d", b, bands[i-1])
			}
			bands = append(bands, int32(b))
		}

		// Before the writer exists, since it needs them at construction.
		contigs, err := collectContigs(cmd, args, vcfRegion)
		if err != nil {
			return err
		}
		if len(contigs) == 0 && vcfVerbose {
			fmt.Fprintf(cmd.ErrOrStderr(),
				"note: no ##contig lines in the input; a VCF exported from this store "+
					"cannot declare its reference\n")
		}

		// Declared before the writer because the writer's rotate hook calls back
		// into it: a shard boundary has to flush the converter's open runs into
		// the shard they belong to, before that shard closes.
		var conv *parquetConverter

		w, err := varstore.NewWriter(vcfToVarstoreOut, varstore.WriterOpts{
			Sink:         sink,
			Codec:        codec,
			RowGroupSize: int64(vcfToVarstoreRowGroupSize),
			Samples:      samples,
			MinDP:        int32(vcfToVarstoreMinDP),
			NoCallable:   vcfToVarstoreNoCallable,
			Program:      buildinfo.String(),
			Command:      buildinfo.CommandLine(),
			Sources:      args,
			Contigs:      contigs,
			Meta:         meta,
			Info:         infoFields,
			Format:       formatFields,
			DepthBands:   bands,
			ShardSites:   int64(vcfToVarstoreShardSites),
			// A run must lie wholly inside one shard, or every locus it covers
			// in an earlier one reads as never assayed rather than reference.
			// The converter already breaks runs at chromosome changes and at
			// depth-band boundaries; a shard boundary is one more of the same.
			BeforeRotate: func() error {
				if conv == nil {
					return nil // no rows written yet, so nothing is open
				}
				return conv.closeRuns()
			},
		})
		if err != nil {
			return err
		}

		started := time.Now()
		conv = &parquetConverter{
			w:          w,
			samples:    samples,
			minDP:      int32(vcfToVarstoreMinDP),
			noCallable: vcfToVarstoreNoCallable,
			runs:       make([]*callableRun, len(samples)),
			bands:      bands,
			format:     formatFields,
			formatBuf:  make(map[string]any, len(formatFields)),
			info:       infoFields,
			infoBuf:    make(map[string]any, len(infoFields)),
			verbose:    vcfVerbose,
			progress:   cmd.ErrOrStderr(),
		}

		for _, path := range args {
			if err := convertOne(cmd, conv, path, samples); err != nil {
				return discarding(w, err)
			}
		}
		if err := conv.finish(); err != nil {
			return discarding(w, err)
		}
		// Before Close, and discarding: this used to run after the store was
		// written, so the failure left all three members on disk -- which then
		// tripped the overwrite guard on the --no-callable retry the message
		// itself asks for.
		if conv.sawDP == 0 && !vcfToVarstoreNoCallable {
			return discarding(w, fmt.Errorf("no DP field found in %s, so callable regions cannot be built\n"+
				"       re-run with --no-callable to accept a store that cannot distinguish\n"+
				"       non-carrier from not-assayed", strings.Join(args, ", ")))
		}
		// Finish closes the members and writes the manifest that marks the store
		// complete; without it the store is unreadable by design. Discard on
		// failure: this was the one error path that returned without cleaning up,
		// and Close is exactly where a full disk shows up -- leaving members that
		// look like a store and then block the retry through the overwrite guard.
		if err := w.Finish(); err != nil {
			return discarding(w, err)
		}

		if vcfVerbose {
			conv.report(cmd.ErrOrStderr(), vcfToVarstoreOut, time.Since(started))
		}
		fmt.Fprintf(cmd.ErrOrStderr(),
			"wrote %s: %d calls, %d sites, %d callable runs over %d samples\n",
			vcfToVarstoreOut, w.NCalls, w.NSites, w.NRegions, len(samples))
		return nil
	},
}

// samplePermutation maps this file's genotype columns onto the canonical
// sample order, and fails if the file does not carry exactly the same people.
//
// Genotype columns are addressed positionally, so getting this wrong does not
// error -- it silently attributes every genotype to the wrong person and
// produces entirely plausible output. Reordering is therefore remapped rather
// than merely tolerated: a bcftools merge or -S reorder is easy to do by
// accident, and remapping turns a silent corruption into a correct result.
func samplePermutation(canonical, got []string, path string) ([]int, bool, error) {
	if len(got) != len(canonical) {
		return nil, false, sampleMismatch(canonical, got, path)
	}
	index := make(map[string]int, len(canonical))
	for i, s := range canonical {
		index[s] = i
	}
	perm := make([]int, len(got))
	reordered := false
	for i, s := range got {
		j, ok := index[s]
		if !ok {
			return nil, false, sampleMismatch(canonical, got, path)
		}
		perm[i] = j
		if i != j {
			reordered = true
		}
	}
	return perm, reordered, nil
}

// sampleMismatch reports which samples differ, since "sample lists differ" on
// a 3,000-sample cohort is not an actionable message.
func sampleMismatch(canonical, got []string, path string) error {
	have := map[string]bool{}
	for _, s := range got {
		have[s] = true
	}
	want := map[string]bool{}
	for _, s := range canonical {
		want[s] = true
	}
	var missing, extra []string
	for _, s := range canonical {
		if !have[s] {
			missing = append(missing, s)
		}
	}
	for _, s := range got {
		if !want[s] {
			extra = append(extra, s)
		}
	}
	msg := fmt.Sprintf("%s does not carry the same samples as the first input (%d vs %d)",
		path, len(got), len(canonical))
	if len(missing) > 0 {
		msg += "\n       missing: " + summariseNames(missing)
	}
	if len(extra) > 0 {
		msg += "\n       unexpected: " + summariseNames(extra)
	}
	return fmt.Errorf("%s", msg)
}

// summariseNames lists a few names rather than thousands.
func summariseNames(names []string) string {
	const max = 6
	if len(names) <= max {
		return strings.Join(names, ", ")
	}
	return fmt.Sprintf("%s and %d more", strings.Join(names[:max], ", "), len(names)-max)
}

// convertOne streams one input into the store.
func convertOne(cmd *cobra.Command, conv *parquetConverter, path string, canonical []string) error {
	src, err := openRecordSource(cmd, path, vcfRegion)
	if err != nil {
		return err
	}
	defer src.close()

	perm, reordered, err := samplePermutation(canonical, src.header.Samples(), path)
	if err != nil {
		return err
	}
	conv.perm = perm

	if conv.verbose {
		note := ""
		if reordered {
			note = "  (sample columns reordered to match the first input)"
		}
		fmt.Fprintf(conv.progress, "reading %s (%d samples)%s\n", path, len(canonical), note)
	}

	for {
		rec, err := src.next()
		if err == io.EOF {
			return nil
		}
		if err != nil {
			return err
		}
		if vcfPassing && rec.IsFiltered() {
			conv.nFiltered++
			continue
		}
		// The store has no way to represent a span of reference, so refuse the
		// moment one turns up rather than silently keeping its first base.
		if rec.IsRefBlock() {
			return gvcfRefBlockError(path, rec)
		}
		if err := conv.checkOrder(rec, path); err != nil {
			return err
		}
		if err := conv.record(rec); err != nil {
			return err
		}
	}
}

// callableRun is an in-progress run of covered sites for one sample.
type callableRun struct {
	start  int32
	last   int32
	nSites int32

	// minDP is the lowest depth seen inside the run, and band is the depth
	// class it belongs to. A run is broken when depth leaves its band, so the
	// bound stays tight: without that, one poorly covered site drags a whole
	// megabase down and every reference call across it inherits the worst
	// moment in the span.
	minDP int32
	band  int
}

// parquetConverter turns a stream of VCF records into store rows, holding only
// one open run per sample so memory does not grow with the input.
type parquetConverter struct {
	w          *varstore.Writer
	samples    []string
	minDP      int32
	noCallable bool
	runs       []*callableRun
	curChrom   string
	sawDP      int64

	// bands are the depth-class boundaries runs are broken at. Empty leaves
	// runs unbanded, which is what stores written before this are.
	bands []int32

	// format are the captured FORMAT fields, and formatBuf is the per-call
	// scratch they are read into. Nil when --format was not given.
	format    []varstore.FormatField
	formatBuf map[string]any

	// info are the captured INFO fields, and infoBuf is the per-record scratch
	// they are read into. Nil when --info was not given.
	info    []varstore.InfoField
	infoBuf map[string]any

	// perm maps the current file's genotype columns onto canonical sample
	// indices; nil means the columns are already in canonical order.
	perm []int

	// ordering cursor, to keep the concatenation coordinate sorted
	lastPos   int32
	seenChrom map[string]bool

	// verbose reporting
	verbose       bool
	progress      io.Writer
	nRecords      int64
	nFiltered     int64
	nMultiAllelic int64
	nExtraRows    int64
	nBlockAlts    int64
	sawGQ         int64
	sawAD         int64
	nGenotypes    int64
	nBelowDP      int64
	nNoCall       int64
	chroms        []string
}

// note records a per-chromosome transition for the verbose summary.
//
// The first coordinate is reported because tick only speaks every 100,000 records,
// so a contig holding fewer than that would otherwise name itself and then go
// silent -- and on a per-chromosome callset that is most of them.
func (c *parquetConverter) note(chrom string, pos int) {
	c.chroms = append(c.chroms, chrom)
	if c.verbose {
		fmt.Fprintf(c.progress, "  %s: starting at %d\n", chrom, pos)
	}
}

// tick emits periodic progress, which matters because a whole-chromosome
// conversion streams for minutes with nothing else to show for it.
//
// The coordinate is the point of it. A record count says the process is alive but
// not how far along it is, and the two are not interchangeable on a sparse input:
// counts advance at the rate variants are called, so a quiet centromere or a
// callset thinned by --passing looks identical to a stall. A position is checkable
// against the contig length, and it is what a run interrupted midway has to be
// resumed or reasoned about from.
func (c *parquetConverter) tick(pos int) {
	const every = 100_000
	if c.verbose && c.nRecords%every == 0 {
		fmt.Fprintf(c.progress, "  %s:%d  %d records, %d calls so far\n",
			c.curChrom, pos, c.nRecords, c.w.NCalls)
	}
}

// record splits one VCF record into per-allele calls, a catalog entry per
// allele, and callable-run bookkeeping.
func (c *parquetConverter) record(rec *vcf.VcfRecord) error {
	chrom := rec.Chrom
	if chrom != c.curChrom {
		// Runs cannot span chromosomes.
		if err := c.closeRuns(); err != nil {
			return err
		}
		c.curChrom = chrom
		c.note(chrom, rec.Pos)
	}
	// BEFORE ANY RUN IS EXTENDED INTO THIS RECORD. A run is extended to this
	// position and only then is the site written, so a rotation discovered
	// while writing the site would emit a run reaching into the shard that is
	// opening -- and every sample would read as never assayed at the first site
	// of every shard. Asking first is what keeps runs inside one shard.
	if c.w.WouldRotate(chrom) {
		if err := c.closeRuns(); err != nil {
			return err
		}
		// ROTATE HERE, not later. Letting the writer discover the boundary while
		// this record's sites are written would put this record's runs in the
		// shard that is closing, and every locus they cover would read as never
		// assayed.
		if err := c.w.Rotate(); err != nil {
			return err
		}
	}

	c.nRecords++
	c.tick(rec.Pos)

	alts := rec.Alt()
	pos := int32(rec.Pos)
	// How far the source record reached. The store's writer would otherwise derive
	// this from len(REF) alone, which is right for a plain variant but understates a
	// symbolic ALT or a gVCF block -- and a region query then misses the record.
	refEnd := int32(rec.RefSpanEnd())
	nAlts := len(alts)
	// A gVCF-derived callset writes the block allele beside a real one
	// ("G,<NON_REF>"). The record is a variant record -- a pure block was refused
	// upstream in convertOne -- but <NON_REF> is not an allele anyone carries, so it
	// is masked out the way any non-focal alternate is. Counted for the -v summary
	// because a store missing it is a question worth being able to answer. Indices
	// into alts stay as the record wrote them, since acCounts is addressed by GT
	// allele number.
	nReal := 0
	for _, alt := range alts {
		if !vcf.IsRefBlockAlt(alt) {
			nReal++
		}
	}
	c.nBlockAlts += int64(nAlts - nReal)
	if nReal > 1 {
		c.nMultiAllelic++
		c.nExtraRows += int64(nReal - 1)
	}
	carriers := make([]int32, nAlts)
	acCounts := make([]int32, nAlts)
	var an int32
	var nLowDP, nCalled int32

	n := rec.NumSamples()
	if n > len(c.samples) {
		n = len(c.samples)
	}
	for i := 0; i < n; i++ {
		sf, err := varstore.ReadSample(rec, i)
		if err != nil {
			return fmt.Errorf("%w (%s:%d)", err, rec.Chrom, rec.Pos)
		}
		c.nGenotypes++
		if sf.DP != varstore.Missing {
			c.sawDP++
		}
		if sf.GQ != varstore.Missing {
			c.sawGQ++
		}
		if sf.AD != "" {
			c.sawAD++
		}
		if !varstore.HasCall(sf.GT) {
			c.nNoCall++
		} else if sf.DP != varstore.Missing && sf.DP < c.minDP {
			c.nBelowDP++
		}

		// Allele counts come straight from GT and are deliberately outside the
		// --no-callable guard below: AC/AN are properties of the genotypes, not
		// of coverage, so they stay meaningful even for a source with no DP.
		an += varstore.AddAlleleCounts(sf.GT, acCounts)

		// Coverage bookkeeping is per site, not per allele. A site counts as
		// callable only when the caller actually made a call there AND depth
		// clears the threshold; "./." at high depth is a declined call, not a
		// covered one.
		if !c.noCallable {
			si := c.sampleAt(i)
			if varstore.HasCall(sf.GT) && sf.DP != varstore.Missing && sf.DP >= c.minDP {
				nCalled++
				band := varstore.DepthBand(c.bands, sf.DP)
				r := c.runs[si]
				// A run spans ONE depth class. Leaving the band closes it and
				// opens another, which is what keeps each run's MinDP a bound
				// on the whole of it rather than on its worst moment.
				if r != nil && band != r.band {
					if err := c.emitRun(si); err != nil {
						return err
					}
					r = nil
				}
				if r != nil {
					r.last = pos
					r.nSites++
					if sf.DP < r.minDP {
						r.minDP = sf.DP
					}
				} else {
					c.runs[si] = &callableRun{
						start: pos, last: pos, nSites: 1, minDP: sf.DP, band: band,
					}
				}
			} else {
				nLowDP++
				if err := c.emitRun(si); err != nil {
					return err
				}
			}
		}

		name := c.samples[c.sampleAt(i)]
		for j, alt := range alts {
			if vcf.IsRefBlockAlt(alt) {
				continue
			}
			call, ok := varstore.CallFor(rec, name, sf, j+1, alt)
			if !ok {
				continue
			}
			call.RefEnd = refEnd
			carriers[j]++
			if len(c.format) == 0 {
				if err := c.w.WriteCall(call); err != nil {
					return err
				}
			} else {
				// j is the ALT index, which is what a Number=A field is indexed
				// by; the sample index is the record's column, not the canonical
				// one, since that is what rec.Sample reads.
				captureFormat(c.formatBuf, rec, i, c.format, j)
				if err := c.w.WriteCallFormat(call, c.formatBuf); err != nil {
					return err
				}
			}
		}
	}

	for j, alt := range alts {
		if vcf.IsRefBlockAlt(alt) {
			continue
		}
		site := varstore.Site{
			Chrom:     chrom,
			Pos:       pos,
			Ref:       rec.Ref,
			Alt:       alt,
			RefEnd:    refEnd,
			AC:        acCounts[j],
			AN:        an,
			NCarriers: carriers[j],
			NLowDP:    nLowDP,
			NCalled:   nCalled,
		}
		if len(c.info) == 0 {
			if err := c.w.WriteSite(site); err != nil {
				return err
			}
			continue
		}
		// j is the ALT index, which is exactly what a Number=A field is indexed
		// by: one record with three ALTs becomes three site rows, and each must
		// take its own value rather than the record's first.
		captureInfo(c.infoBuf, rec, c.info, j)
		if err := c.w.WriteSiteInfo(site, c.infoBuf); err != nil {
			return err
		}
	}
	return nil
}

// emitRun writes and clears sample i's open run, if any.
func (c *parquetConverter) emitRun(i int) error {
	r := c.runs[i]
	if r == nil {
		return nil
	}
	c.runs[i] = nil
	return c.w.WriteRegion(varstore.CalledSiteRun{
		SampleID: c.samples[i],
		Chrom:    c.curChrom,
		Start:    r.start,
		End:      r.last,
		NSites:   r.nSites,
		MinDP:    r.minDP,
	})
}

// closeRuns flushes every open run, at a chromosome change or end of input.
func (c *parquetConverter) closeRuns() error {
	for i := range c.runs {
		if err := c.emitRun(i); err != nil {
			return err
		}
	}
	return nil
}

// finish flushes any runs still open at end of input.
func (c *parquetConverter) finish() error { return c.closeRuns() }

// report writes the verbose conversion summary.
//
// The field-presence section is the part worth having. A gate can only act on a
// field the data carries, and --min-gq over GQ-less input admits everything
// rather than rejecting it. That is deliberate -- absent quality is not evidence
// of poor quality -- but it means a filter can silently do nothing, so a store
// built from such input should say so at the point it is created.
func (c *parquetConverter) report(out io.Writer, base string, elapsed time.Duration) {
	pct := func(n int64) string {
		if c.nGenotypes == 0 {
			return "0%"
		}
		return fmt.Sprintf("%.1f%%", 100*float64(n)/float64(c.nGenotypes))
	}

	fmt.Fprintf(out, "\ninput\n")
	fmt.Fprintf(out, "  records read          %d\n", c.nRecords)
	if c.nFiltered > 0 {
		fmt.Fprintf(out, "  skipped (--passing)   %d\n", c.nFiltered)
	}
	fmt.Fprintf(out, "  multiallelic records  %d (split into %d extra rows)\n",
		c.nMultiAllelic, c.nExtraRows)
	if c.nBlockAlts > 0 {
		fmt.Fprintf(out, "  <NON_REF> alts masked %d  (block alleles beside a real one; not stored as alleles)\n",
			c.nBlockAlts)
	}
	fmt.Fprintf(out, "  chromosomes           %s\n", strings.Join(c.chroms, ", "))
	fmt.Fprintf(out, "  genotypes examined    %d\n", c.nGenotypes)

	fmt.Fprintf(out, "\nfields present (a gate can only act on a field the data has)\n")
	for _, f := range []struct {
		name string
		n    int64
		note string
	}{
		{"DP", c.sawDP, "callable runs and --min-dp depend on this"},
		{"GQ", c.sawGQ, "--min-gq depends on this"},
		{"AD", c.sawAD, "per-allele depths"},
	} {
		state := pct(f.n)
		if f.n == 0 {
			state = "ABSENT -- " + f.note + " will have no effect"
		}
		fmt.Fprintf(out, "  %-3s %s\n", f.name, state)
	}

	fmt.Fprintf(out, "\ncoverage at --min-dp %d\n", c.minDP)
	if c.noCallable {
		fmt.Fprintf(out, "  not tracked (--no-callable)\n")
	} else {
		fmt.Fprintf(out, "  no call made          %d  (%s)\n", c.nNoCall, pct(c.nNoCall))
		fmt.Fprintf(out, "  called but below DP   %d  (%s)\n", c.nBelowDP, pct(c.nBelowDP))
	}

	fmt.Fprintf(out, "\noutput\n")
	for _, f := range []struct {
		path string
		n    int64
	}{
		{varstore.CallsPath(base), c.w.NCalls},
		{varstore.SitesPath(base), c.w.NSites},
		{varstore.RegionsPath(base), c.w.NRegions},
	} {
		size := int64(-1)
		if st, err := os.Stat(f.path); err == nil {
			size = st.Size()
		}
		fmt.Fprintf(out, "  %-24s %9d rows  %10d bytes\n", filepath.Base(f.path), f.n, size)
	}
	fmt.Fprintf(out, "  elapsed               %s\n", elapsed.Round(time.Millisecond))
}

func init() {
	f := vcfToVarstoreCmd.Flags()
	f.StringVar(&vcfToVarstoreOut, "out", "", "Store directory, created if needed (DIR/calls.parquet etc)")
	addRegionFlag(vcfToVarstoreCmd)
	f.IntVar(&vcfToVarstoreMinDP, "min-dp", 10, "Minimum DP for a site to count as callable for a sample")
	f.IntVar(&vcfToVarstoreShardSites, "shard-sites", 0,
		"Split each member every N sites, so a locus query reads one small file instead of pruning a large one (0 writes one file per member)")
	f.IntSliceVar(&vcfToVarstoreBands, "depth-bands", []int{10, 20, 50},
		"Depth boundaries at which a callable run is broken, so each run's min_dp bounds the whole of it (empty leaves runs unbanded)")
	f.StringSliceVar(&vcfToVarstoreFormat, "format", nil,
		"FORMAT field to capture onto each ALT call, as its own column (repeatable, comma separated, globs allowed)")
	f.StringSliceVar(&vcfToVarstoreInfo, "info", nil,
		"INFO field to capture into the sites catalog, as its own column (repeatable, comma separated, globs allowed)")
	f.BoolVar(&vcfToVarstoreNoCallable, "no-callable", false, "Accept a source with no DP field; callable regions will be empty")
	addPassingFlag(vcfToVarstoreCmd, "Only convert passing variants")
	f.StringVar(&vcfToVarstoreCompression, "compression", "zstd", "Parquet compression: zstd, snappy, or none")
	f.IntVar(&vcfToVarstoreRowGroupSize, "row-group-size", 250000, "Rows per parquet row group")
	addVerboseFlag(vcfToVarstoreCmd, "Report progress and a conversion summary on stderr")
	f.BoolVar(&vcfToVarstoreForce, "force", false, "Overwrite an existing store at --out")

	// Generated from the library's own list, so the flag names and the keys they
	// write cannot drift apart, and a key added upstream shows up here for free.
	for _, key := range varstore.ReservedMetaKeys {
		p := new(string)
		vcfToVarstoreMetaNamed[key] = p
		f.StringVar(p, "meta-"+key, "", metaFlagHelp[key])
	}
	f.StringArrayVar(&vcfToVarstoreMetaPairs, "meta", nil,
		"Record KEY=VALUE in the store's metadata (repeatable; lowercase key of [a-z0-9_-])")
}

// metaFlagHelp is the one-line help for each reserved metadata key. Keyed by the
// key rather than positional so a reordering upstream cannot mislabel a flag.
var metaFlagHelp = map[string]string{
	varstore.MetaKeyDataset:     "Name of the dataset or release this store was built from",
	varstore.MetaKeyReference:   "Reference assembly the calls were made against (e.g. GRCh38)",
	varstore.MetaKeyCaller:      "Variant caller and version that produced the input (e.g. 'GATK 4.2.6.1')",
	varstore.MetaKeyAccession:   "Study or dataset accession (e.g. phs000000)",
	varstore.MetaKeyURL:         "Where the source data was retrieved from",
	varstore.MetaKeyVersion:     "Release version of the dataset itself",
	varstore.MetaKeyDescription: "Free-text description of the dataset",
}

// sampleAt maps a genotype column of the file being read to its canonical
// sample index.
func (c *parquetConverter) sampleAt(col int) int {
	if c.perm == nil || col >= len(c.perm) {
		return col
	}
	return c.perm[col]
}

// checkOrder enforces that the inputs, concatenated, stay coordinate sorted.
//
// This is one rule serving three purposes: it keeps parquet's per-row-group
// min/max on pos tight, which is what makes locus lookups prune; it catches
// inputs supplied in the wrong order; and it rejects overlapping inputs, which
// would otherwise write duplicate sites and split AC/AN across two rows for the
// same variant.
func (c *parquetConverter) checkOrder(rec *vcf.VcfRecord, path string) error {
	if c.seenChrom == nil {
		c.seenChrom = map[string]bool{}
	}
	chrom, pos := rec.Chrom, int32(rec.Pos)
	if chrom == c.curChrom {
		if pos < c.lastPos {
			return fmt.Errorf("%s is not coordinate sorted at %s:%d (previous record was %d)\n"+
				"       inputs must be sorted, and must not overlap each other",
				path, chrom, pos, c.lastPos)
		}
		c.lastPos = pos
		return nil
	}
	if c.seenChrom[chrom] {
		return fmt.Errorf("%s returns to %s after another chromosome\n"+
			"       inputs must be in coordinate order and must not overlap; "+
			"pass them one chromosome at a time, in order", path, chrom)
	}
	c.seenChrom[chrom] = true
	c.lastPos = pos
	return nil
}

// collectContigs gathers the ##contig lines of every input, so the store can be
// exported back to VCF with a header that says which reference it came from.
//
// The UNION, not the first input's. A whole-genome callset usually ships one VCF
// per chromosome, and such a file often declares only its own contig -- taking the
// first input's list would silently lose every other chromosome, which is exactly
// the case multi-input conversion exists to serve.
//
// Lines are kept verbatim, so lengths and any extra fields survive rather than
// being reconstructed approximately.
func collectContigs(cmd *cobra.Command, inputs []string, region string) ([]string, error) {
	type origin struct {
		line string
		from string
	}
	byID := map[string]origin{}
	var order []string

	for _, path := range inputs {
		src, err := openRecordSource(cmd, path, region)
		if err != nil {
			return nil, err
		}
		names := src.header.ContigNames()
		for _, id := range names {
			def, ok := src.header.ContigDef(id)
			if !ok {
				continue
			}
			line := def.String()
			prev, seen := byID[id]
			if !seen {
				byID[id] = origin{line: line, from: path}
				order = append(order, id)
				continue
			}
			// Same contig declared differently by two inputs. A differing length
			// means the inputs were called against different references, which would
			// make one store out of two incompatible callsets -- refused for the same
			// reason a differing sample set is.
			if prev.line != line {
				src.close()
				return nil, fmt.Errorf("inputs disagree about contig %s:\n"+
					"       %s: %s\n"+
					"       %s: %s\n"+
					"       a differing length means these were called against different "+
					"references", id, prev.from, prev.line, path, line)
			}
		}
		src.close()
	}

	out := make([]string, 0, len(order))
	for _, id := range order {
		out = append(out, byID[id].line)
	}
	return out, nil
}

// discarding removes the half-written store and returns why the conversion
// failed, mentioning the cleanup only when that failed too.
//
// A failed removal matters enough to say: what survives looks like a store, so
// the overwrite guard will refuse to replace it, and the retry the error asks
// for cannot run until it is deleted by hand.
func discarding(w *varstore.Writer, err error) error {
	if derr := w.Discard(); derr != nil {
		return fmt.Errorf("%w\n       the partial store could not be removed: %v", err, derr)
	}
	return err
}

// writableSchemes names the remote destinations this build can reach, for the
// error a bare "cannot write there" would otherwise leave unexplained.
//
// A transport is linked in by import, so which ones exist is a property of the
// binary rather than of the arguments -- and "s3:// is not supported" and
// "s3:// support was not compiled into this build" are different problems with
// different fixes.
func writableSchemes() string {
	schemes := varstore.SinkSchemes()
	if len(schemes) == 0 {
		return ", and no remote transports are linked into this build"
	}
	for i, s := range schemes {
		schemes[i] = s + "://"
	}
	return ", or " + strings.Join(schemes, ", ")
}
