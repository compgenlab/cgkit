package samcmd

import (
	"fmt"
	"os"
	"strings"

	"github.com/compgenlab/hts/htsio"
	_ "github.com/compgenlab/hts/htsio/bam"
	_ "github.com/compgenlab/hts/htsio/cram"
	_ "github.com/compgenlab/hts/htsio/sam"
	"github.com/compgenlab/hts/seqio"
	"github.com/compgenlab/hts/support/sequtils"
	"github.com/compgenlab/hts/support/stringutils"
	"github.com/spf13/cobra"
)

var samToFastaCmd = &cobra.Command{
	GroupID: "samcmd",
	Use:     "sam-tofasta <input.bam> [output]",
	Short:   "Convert SAM/BAM/CRAM reads to FASTA",
	Long:    "Write SAM/BAM/CRAM reads to FASTA. Output file is optional; defaults to stdout.",
	RunE: func(cmd *cobra.Command, args []string) error {
		return runSamToFastx(cmd, args, false, samToFastaReaderFlags, samToFastaWriteTags)
	},
}

var samToFastqCmd = &cobra.Command{
	GroupID: "samcmd",
	Use:     "sam-tofastq <input.bam> [output]",
	Short:   "Convert SAM/BAM/CRAM reads to FASTQ",
	Long:    "Write SAM/BAM/CRAM reads to FASTQ. Output file is optional; defaults to stdout.",
	RunE: func(cmd *cobra.Command, args []string) error {
		return runSamToFastx(cmd, args, true, samToFastqReaderFlags, samToFastqWriteTags)
	},
}

// runSamToFastx converts SAM/BAM/CRAM reads to FASTA (fastq=false) or FASTQ
// (fastq=true). Reverse-strand reads are reverse-complemented. Tags listed in
// writeTags are emitted into the record comment as tab-delimited SAM
// tag:type:value fields; tags absent from a record are omitted.
func runSamToFastx(cmd *cobra.Command, args []string, fastq bool, readerFlags samReaderFlags, writeTags []string) error {
	if len(args) == 0 {
		cmd.Help()
		return nil
	}

	opts, err := readerFlags.buildReaderOpts()
	if err != nil {
		return err
	}

	tags := normalizeTags(writeTags)

	inputFile := args[0]
	reader, err := htsio.NewSamReader(inputFile, opts)
	if err != nil {
		return err
	}
	defer reader.Close()

	var fastaWriter *seqio.FastaWriter
	var fastqWriter *seqio.FastqWriter

	if len(args) > 1 {
		outputFile := args[1]
		if fastq {
			fastqWriter, err = seqio.OpenFastqWriter(outputFile)
			if err != nil {
				return fmt.Errorf("open FASTQ output: %w", err)
			}
		} else {
			fastaWriter, err = seqio.OpenFastaWriter(outputFile)
			if err != nil {
				return fmt.Errorf("open FASTA output: %w", err)
			}
		}
	} else {
		if fastq {
			fastqWriter = seqio.NewFastqWriter(os.Stdout)
		} else {
			fastaWriter = seqio.NewFastaWriter(os.Stdout)
		}
	}
	if fastqWriter != nil {
		defer fastqWriter.Close()
	}
	if fastaWriter != nil {
		defer fastaWriter.Close()
	}

	for rec, err := range reader.Records() {
		if err != nil {
			return err
		}

		seq := rec.Seq
		qual := rec.Qual

		// Reverse complement if on reverse strand.
		if rec.IsReverse() {
			seq = sequtils.ReverseComplement(seq)
			qual = stringutils.ReverseString(qual)
		}

		comment := buildTagComment(rec, tags)

		if fastq {
			if err := fastqWriter.WriteRecord(rec.ReadName, comment, seq, qual); err != nil {
				return fmt.Errorf("write FASTQ: %w", err)
			}
		} else {
			if err := fastaWriter.WriteRecord(rec.ReadName, comment, seq); err != nil {
				return fmt.Errorf("write FASTA: %w", err)
			}
		}
	}

	return nil
}

// normalizeTags expands a repeatable, comma-separated --write-tag flag into an
// ordered list of two-letter tag codes, preserving the order given and dropping
// empty entries.
func normalizeTags(writeTags []string) []string {
	var tags []string
	for _, v := range writeTags {
		for _, t := range strings.Split(v, ",") {
			if t = strings.TrimSpace(t); t != "" {
				tags = append(tags, t)
			}
		}
	}
	return tags
}

// buildTagComment formats the requested tags present on rec as tab-delimited
// SAM tag:type:value fields, matching SamRecord.String(). Tags absent from the
// record are omitted; the result is empty if none are present.
func buildTagComment(rec *htsio.SamRecord, tags []string) string {
	if len(tags) == 0 {
		return ""
	}
	var parts []string
	for _, tag := range tags {
		if t, ok := rec.Tags[tag]; ok {
			parts = append(parts, fmt.Sprintf("%s:%c:%s", tag, t.Type, t.Value))
		}
	}
	return strings.Join(parts, "\t")
}

var (
	samToFastaReaderFlags samReaderFlags
	samToFastaWriteTags   []string
	samToFastqReaderFlags samReaderFlags
	samToFastqWriteTags   []string
)

func init() {
	samToFastaCmd.Flags().StringArrayVar(&samToFastaWriteTags, "write-tag", nil, "SAM tag to write into the read comment (repeatable; comma-separated)")
	samToFastaReaderFlags.register(samToFastaCmd)

	samToFastqCmd.Flags().StringArrayVar(&samToFastqWriteTags, "write-tag", nil, "SAM tag to write into the read comment (repeatable; comma-separated)")
	samToFastqReaderFlags.register(samToFastqCmd)
}
