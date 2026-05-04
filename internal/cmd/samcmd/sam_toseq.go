package samcmd

import (
	"fmt"
	"os"

	"github.com/compgen-io/cgltk/htsio"
	"github.com/compgen-io/cgltk/seqio"
	"github.com/compgen-io/cgltk/support/sequtils"
	"github.com/compgen-io/cgltk/support/stringutils"
	"github.com/spf13/cobra"
)

var samToSeqCmd = &cobra.Command{
	GroupID: "samcmd",
	Use:     "sam-toseq <input.bam> [output]",
	Short:   "Convert SAM/BAM/CRAM reads to FASTA or FASTQ",
	Long:    "Write SAM/BAM/CRAM reads to FASTA or FASTQ. Output file is optional; defaults to stdout.",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}

		if !samToSeqFasta && !samToSeqFastq {
			return fmt.Errorf("at least one of --fasta or --fastq must be specified")
		}

		opts, err := samToSeqReaderFlags.buildReaderOpts()
		if err != nil {
			return err
		}

		inputFile := args[0]
		reader, err := htsio.NewSamReader(inputFile, opts)
		if err != nil {
			return err
		}
		defer reader.Close()

		// Open output writers.
		var fastaWriter *seqio.FastaWriter
		var fastqWriter *seqio.FastqWriter

		if len(args) > 1 {
			outputFile := args[1]
			if samToSeqFasta {
				fastaWriter, err = seqio.OpenFastaWriter(outputFile)
				if err != nil {
					return fmt.Errorf("open FASTA output: %w", err)
				}
				defer fastaWriter.Close()
			} else {
				fastqWriter, err = seqio.OpenFastqWriter(outputFile)
				if err != nil {
					return fmt.Errorf("open FASTQ output: %w", err)
				}
				defer fastqWriter.Close()
			}
		} else {
			if samToSeqFasta {
				fastaWriter = seqio.NewFastaWriter(os.Stdout)
				defer fastaWriter.Close()
			} else {
				fastqWriter = seqio.NewFastqWriter(os.Stdout)
				defer fastqWriter.Close()
			}
		}

		for rec, err := range reader.Records() {
			if err != nil {
				return err
			}

			seq := rec.Seq
			qual := rec.Qual

			// Reverse complement if on reverse strand.
			if rec.IsReverse() {
				seq = sequtils.ReverseCompliment(seq)
				qual = stringutils.ReverseString(qual)
			}

			if fastaWriter != nil {
				if err := fastaWriter.WriteRecord(rec.ReadName, "", seq); err != nil {
					return fmt.Errorf("write FASTA: %w", err)
				}
			}
			if fastqWriter != nil {
				if err := fastqWriter.WriteRecord(rec.ReadName, "", seq, qual); err != nil {
					return fmt.Errorf("write FASTQ: %w", err)
				}
			}
		}

		return nil
	},
}

var (
	samToSeqFasta       bool
	samToSeqFastq       bool
	samToSeqReaderFlags samReaderFlags
)

func init() {
	samToSeqCmd.Flags().BoolVar(&samToSeqFasta, "fasta", false, "Output in FASTA format")
	samToSeqCmd.Flags().BoolVar(&samToSeqFastq, "fastq", false, "Output in FASTQ format")
	samToSeqReaderFlags.register(samToSeqCmd)
}
