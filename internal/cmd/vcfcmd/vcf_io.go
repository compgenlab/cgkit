package vcfcmd

import (
	"fmt"
	"io"
	"iter"
	"math"
	"os"
	"strings"

	"github.com/compgenlab/cghts/htsio"
	"github.com/compgenlab/cghts/vcf"
	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/compgenlab/cgkit/internal/cmdio"
	"github.com/spf13/cobra"
)

const sinceVersion = "v0.4.0"

// stampVcfProvenance updates the header's ##fileDate and appends cgkit command
// and version provenance lines, recording how the output file was produced. It
// is called by any command that writes a VCF file.
func stampVcfProvenance(h *vcf.VcfHeader, cmdName string) {
	h.SetFileDate(buildinfo.Date())
	h.AddLine("##cgkit_" + cmdName + "Command=" + buildinfo.CommandLine())
	h.AddLine("##cgkit_" + cmdName + "Version=" + buildinfo.String())
}

// openVcfInput opens a streaming VCF reader for filename, reading from stdin
// when filename is "-". Input is transparently gunzipped, and may be a local
// path or a remote locator (http(s)://, s3://).
//
// Note that a streaming read of a remote object transfers the whole thing;
// there is no index to skip with. Commands taking --region seek instead.
func openVcfInput(cmd *cobra.Command, filename string) (*vcf.VcfReader, error) {
	if filename == "-" {
		return vcf.NewVcfReader(cmd.InOrStdin())
	}
	r, err := vcf.OpenVcfFile(cmd.Context(), filename)
	if err != nil {
		// The locator is worth repeating: for a remote read the underlying
		// error is an HTTP status or an SDK message with no path in it.
		return nil, fmt.Errorf("opening %s: %w", filename, err)
	}
	return r, nil
}

// openOutput returns the writer for output, using stdout when output is "" or
// "-". The returned closer is nil when writing to stdout.
func openOutput(cmd *cobra.Command, output string) (io.Writer, func() error, error) {
	out, err := cmdio.CreateOutput(cmd, "-o/--output", output)
	if err != nil {
		return nil, nil, err
	}
	if output == "" || output == "-" {
		return out.W, nil, nil
	}
	return out.W, out.Close, nil
}

// recordSource is a uniform record stream over either a plain (streaming) VCF or
// a tabix-indexed region query. next returns io.EOF when exhausted.
type recordSource struct {
	header *vcf.VcfHeader
	next   func() (*vcf.VcfRecord, error)
	close  func() error
}

// openRecordSource opens filename for reading. When region is non-empty it must
// name a tabix-indexed file (not stdin); records are limited to that region
// (1-based inclusive "chrom:start-end", or bare "chrom"). Otherwise the file is
// read as a stream.
func openRecordSource(cmd *cobra.Command, filename, region string) (*recordSource, error) {
	if region != "" {
		if filename == "-" {
			return nil, fmt.Errorf("--region requires an indexed VCF file, not stdin")
		}
		ref, start, end, err := htsio.ParseRegion(region)
		if err != nil {
			return nil, err
		}
		if end < 0 {
			end = math.MaxInt32
		}
		ir, err := vcf.OpenIndexedVcfReader(cmd.Context(), filename)
		if err != nil {
			return nil, err
		}
		header, err := ir.Header()
		if err != nil {
			ir.Close()
			return nil, err
		}
		seq, err := ir.Query(ref, start, end)
		if err != nil {
			ir.Close()
			return nil, err
		}
		next, stop := iter.Pull2(seq)
		return &recordSource{
			header: header,
			next: func() (*vcf.VcfRecord, error) {
				rec, qerr, ok := next()
				if !ok {
					return nil, io.EOF
				}
				if qerr != nil {
					return nil, qerr
				}
				return rec, nil
			},
			close: func() error {
				stop()
				return ir.Close()
			},
		}, nil
	}

	reader, err := openVcfInput(cmd, filename)
	if err != nil {
		return nil, err
	}
	header, err := reader.Header()
	if err != nil {
		reader.Close()
		return nil, err
	}
	return &recordSource{
		header: header,
		next:   reader.NextRecord,
		close: func() error {
			reader.Close()
			return nil
		},
	}, nil
}

// plural returns "s" when n is not 1, for error messages that name a list.
func plural(n int) string {
	if n == 1 {
		return ""
	}
	return "s"
}

// resolveSampleIndex looks up a sample by name or by 1-based position.
//
// The bounds check is the point. header.SampleIndex resolves a name it does not
// know to -1, but a *numeric* name to n-1 with no upper bound at all, so
// "--sample 9" against a 3-sample file returned 8 and passed any `idx < 0`
// guard. What followed was a per-record "sample index 8 out of range" from deep
// inside the reader, once per variant, rather than one error naming the flag.
func resolveSampleIndex(header *vcf.VcfHeader, flag, name string) (int, error) {
	samples := header.Samples()
	idx := header.SampleIndex(name)
	if idx < 0 || idx >= len(samples) {
		return 0, fmt.Errorf("%s: no such sample %q\n  this file has: %s",
			flag, name, strings.Join(samples, ", "))
	}
	return idx, nil
}

// vcfStream describes a single-input, single-output VCF transform: read every
// record, optionally change or drop it, write the rest.
type vcfStream struct {
	// name is the command name, for the provenance stamp.
	name string
	// in and out are the input locator and the -o value.
	in, out string
	// header, when set, adjusts the header before provenance is stamped onto it.
	header func(*vcf.VcfHeader) error
	// record returns false to drop the record. It may modify it in place.
	record func(*vcf.VcfRecord) (bool, error)
	// write, when set, replaces WriteRecord -- vcf-reorder rewrites the sample
	// columns as a raw line rather than through the record.
	write func(*vcf.VcfWriter, *vcf.VcfRecord) error
}

// runVcfStream runs a vcfStream to completion.
//
// This was eleven near-identical copies of open-header-stamp-write-loop-close,
// and they shared a bug that is the real reason to have one of them. The writer
// was never deferred and the closing tail sat *after* the record loop, so every
// error return inside the loop skipped it: the file descriptor leaked, and for a
// bgzip output the file was left with no BGZF EOF block. That last part is the
// dangerous half -- a truncated BGZF is detectably broken, but only if something
// checks, and in the meantime it is a file sitting where a result belongs.
//
// So a failure closes the writer and removes what it wrote. vcf-toparquet
// already did this (see discarding); everything else did not.
func runVcfStream(cmd *cobra.Command, s vcfStream) (err error) {
	reader, err := openVcfInput(cmd, s.in)
	if err != nil {
		return err
	}
	defer reader.Close()

	header, err := reader.Header()
	if err != nil {
		return err
	}
	if s.header != nil {
		if err := s.header(header); err != nil {
			return err
		}
	}
	stampVcfProvenance(header, s.name)

	writer, closeFn, err := openVcfWriter(cmd, s.out)
	if err != nil {
		return err
	}
	// Set once the stream is fully written and only the closer remains. A
	// failure before that point left a partial file and is removed; a failure
	// from the closer is not, because --tbi refuses to index unsorted output
	// *after* the VCF is complete, and that VCF is valid -- the error names it
	// and says how to sort it. Discarding it would destroy a good result over a
	// missing index.
	closing := false
	defer func() {
		if err == nil || closing {
			return
		}
		writer.Close()
		discardPartialVcf(s.out)
	}()

	if err := writer.WriteHeader(header); err != nil {
		return err
	}
	for {
		rec, err := reader.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		keep, err := s.record(rec)
		if err != nil {
			return err
		}
		if !keep {
			continue
		}
		write := writer.WriteRecord
		if s.write != nil {
			write = func(r *vcf.VcfRecord) error { return s.write(writer, r) }
		}
		if err := write(rec); err != nil {
			return err
		}
	}
	closing = true
	if closeFn != nil {
		return closeFn()
	}
	return writer.Close()
}

// discardPartialVcf removes the half-written output of a failed run, and its
// index if one was requested.
//
// Only a file this process created and truncated is touched -- stdout has
// nothing to remove, and the removal errors are ignored because the failure that
// got us here is the one worth reporting.
func discardPartialVcf(out string) {
	if out == "" || out == "-" {
		return
	}
	os.Remove(out)
	if vcfTbi {
		os.Remove(out + ".tbi")
	}
}

// Flags shared by many commands, backed by one variable each.
//
// This follows vcfTbi: only one subcommand runs per process, so a variable per
// command bought nothing and cost something. --passing had fourteen of them and
// --region five, each registered with its own hand-written help string, and
// every one had to be added by hand to the reset list in vcf_commands_test.go
// or it leaked into the next test.
var (
	// vcfPassing backs --passing. The phrasing differs per command -- "output",
	// "export", "convert" -- so the help text stays a caller's argument even
	// though the flag does not.
	vcfPassing bool

	// vcfRegion backs --region. All five declarations were character-identical.
	vcfRegion string

	// vcfVerbose backs -v/--verbose. What gets reported differs per command, so
	// the help text is still the caller's.
	vcfVerbose bool
)

// addVerboseFlag registers -v/--verbose with command-specific help.
func addVerboseFlag(cmd *cobra.Command, help string) {
	cmd.Flags().BoolVarP(&vcfVerbose, "verbose", "v", false, help)
}

// addOutputFlag registers -o/--output for a command writing a tabular stream
// rather than a VCF. The VCF writers use addVcfOutputFlags, whose help has to
// mention bgzip and --tbi; these had six copies of one string and two lone
// variants, which is the kind of drift a shared registrar exists to stop.
func addOutputFlag(cmd *cobra.Command, out *string) {
	cmd.Flags().StringVarP(out, "output", "o", "-", "Output filename (- for stdout)")
}

// addPassingFlag registers --passing with command-specific help.
func addPassingFlag(cmd *cobra.Command, help string) {
	cmd.Flags().BoolVar(&vcfPassing, "passing", false, help)
}

// addRegionFlag registers --region, which needs a tabix-indexed input.
func addRegionFlag(cmd *cobra.Command) {
	cmd.Flags().StringVar(&vcfRegion, "region", "",
		"Only variants in this 1-based region (chrom:start-end, or chrom); requires a tabix-indexed file")
}
