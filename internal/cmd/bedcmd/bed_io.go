package bedcmd

import (
	"io"

	"github.com/compgenlab/cghts/bed"
	"github.com/spf13/cobra"
)

// openBedInput opens a BED reader for filename, reading from stdin when
// filename is "-". Input is transparently gunzipped.
func openBedInput(cmd *cobra.Command, filename string) (*bed.BedReader, error) {
	if filename == "-" {
		return bed.NewBedReader(cmd.InOrStdin())
	}
	return bed.NewBedFile(filename)
}

// openBedOutput opens a BED writer for output, writing to stdout when output is
// "" or "-". A filename ending in ".gz" is gzip-compressed.
func openBedOutput(cmd *cobra.Command, output string, opts *bed.BedWriterOpts) (*bed.BedWriter, error) {
	if output == "" || output == "-" {
		return bed.NewBedWriter(cmd.OutOrStdout(), opts), nil
	}
	return bed.OpenBedWriter(output, opts)
}

// streamBed reads every record from the input file and writes it through the
// writer configured by opts. It is the common body of the simple
// reader-to-writer BED commands.
func streamBed(cmd *cobra.Command, filename string, opts *bed.BedWriterOpts, output string) error {
	reader, err := openBedInput(cmd, filename)
	if err != nil {
		return err
	}
	defer reader.Close()

	writer, err := openBedOutput(cmd, output, opts)
	if err != nil {
		return err
	}

	for {
		rec, err := reader.NextRecord()
		if err != nil {
			if err == io.EOF {
				break
			}
			return err
		}
		if err := writer.WriteRecord(rec); err != nil {
			return err
		}
	}
	return writer.Close()
}
