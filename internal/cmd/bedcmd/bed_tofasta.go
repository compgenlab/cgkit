package bedcmd

import (
	"fmt"
	"io"

	"github.com/compgenlab/cghts/bed"
	"github.com/compgenlab/cghts/seqio"
	"github.com/compgenlab/cghts/support/sequtils"
	"github.com/spf13/cobra"
)

var (
	bedToFastaOutput       string
	bedToFastaIgnoreStrand bool
	bedToFastaWrap         int
)

var bedToFastaCmd = &cobra.Command{
	GroupID:     "bedcmd",
	Annotations: map[string]string{"since": "v0.3.1"},
	Use:         "bed-tofasta <input.bed> <ref.fa>",
	Short:       "Extract FASTA sequences based on BED coordinates",
	Long: `Extract FASTA sequences based on BED coordinates.

The reference FASTA must be indexed (samtools faidx; the .fai may sit alongside
a bgzip-compressed FASTA with a .gzi index).

By default minus-strand regions are reverse-complemented and their coordinates
are reported high-to-low in the sequence name. With --ignore-strand, the strand
is ignored and the plus-strand sequence is always returned.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) < 2 {
			cmd.Help()
			return nil
		}

		reader, err := openBedInput(cmd, args[0])
		if err != nil {
			return err
		}
		defer reader.Close()

		ref, err := seqio.NewIndexedFastaReader(args[1])
		if err != nil {
			return err
		}
		defer ref.Close()

		opts := seqio.NewFastaWriterOpts().Wrap(bedToFastaWrap)
		var writer *seqio.FastaWriter
		if bedToFastaOutput == "" || bedToFastaOutput == "-" {
			writer = seqio.NewFastaWriter(cmd.OutOrStdout(), opts)
		} else {
			writer, err = seqio.OpenFastaWriter(bedToFastaOutput, opts)
			if err != nil {
				return err
			}
		}

		for {
			rec, err := reader.NextRecord()
			if err != nil {
				if err == io.EOF {
					break
				}
				return err
			}

			bases, err := ref.GetSequenceRange(rec.Ref, rec.Start, rec.End)
			if err != nil {
				return err
			}
			seq := string(bases)

			var name string
			if !bedToFastaIgnoreStrand && rec.Strand == bed.StrandMinus {
				name = fmt.Sprintf("%s|%s:%d-%d", rec.Name, rec.Ref, rec.End, rec.Start)
				seq = sequtils.ReverseComplement(seq)
			} else {
				name = fmt.Sprintf("%s|%s:%d-%d", rec.Name, rec.Ref, rec.Start, rec.End)
			}

			if err := writer.WriteRecord(name, "", seq); err != nil {
				return err
			}
		}
		return writer.Close()
	},
}

func init() {
	bedToFastaCmd.Flags().StringVarP(&bedToFastaOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
	bedToFastaCmd.Flags().BoolVar(&bedToFastaIgnoreStrand, "ignore-strand", false, "Ignore strand (always return the + strand sequence)")
	bedToFastaCmd.Flags().IntVar(&bedToFastaWrap, "wrap", 0, "Wrap output sequence at N bases (0 = no wrapping)")
}
