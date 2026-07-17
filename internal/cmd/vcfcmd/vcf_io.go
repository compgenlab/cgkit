package vcfcmd

import (
	"fmt"
	"io"
	"iter"
	"math"
	"os"

	"github.com/compgenlab/cghts/htsio"
	"github.com/compgenlab/cghts/vcf"
	"github.com/compgenlab/cgkit/internal/buildinfo"
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
// when filename is "-". Input is transparently gunzipped.
func openVcfInput(cmd *cobra.Command, filename string) (*vcf.VcfReader, error) {
	if filename == "-" {
		return vcf.NewVcfReader(cmd.InOrStdin())
	}
	return vcf.NewVcfFile(filename)
}

// openOutput returns the writer for output, using stdout when output is "" or
// "-". The returned closer is nil when writing to stdout.
func openOutput(cmd *cobra.Command, output string) (io.Writer, func() error, error) {
	if output == "" || output == "-" {
		return cmd.OutOrStdout(), nil, nil
	}
	f, err := os.Create(output)
	if err != nil {
		return nil, nil, err
	}
	return f, f.Close, nil
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
		ir, err := vcf.NewIndexedVcfReader(filename)
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
