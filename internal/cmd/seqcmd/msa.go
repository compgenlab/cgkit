package seqcmd

import (
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/compgen-io/cgltk/align"
	"github.com/compgen-io/cgltk/seqio"
	"github.com/spf13/cobra"
)

var msaCmd = &cobra.Command{
	GroupID: "seqcmd",
	Use:     "seq-consensus-msa <input.fasta|fastq>",
	Short:   "Multiple sequence alignment via incremental consensus",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}

		// Detect file format by extension
		filename := args[0]
		isFastq := strings.HasSuffix(filename, ".fq") ||
			strings.HasSuffix(filename, ".fastq") ||
			strings.HasSuffix(filename, ".fq.gz") ||
			strings.HasSuffix(filename, ".fastq.gz")

		var reader seqio.SeqReader
		var err error
		if isFastq {
			reader, err = seqio.NewFastqFile(filename)
		} else {
			reader, err = seqio.NewFastaFile(filename)
		}
		if err != nil {
			return fmt.Errorf("opening input: %w", err)
		}

		// Load all sequences into memory
		var seqs []seqio.SeqQual
		var names []string
		for {
			rec, err := reader.NextSeq()
			if err != nil {
				if err == io.EOF {
					break
				}
				return fmt.Errorf("reading input: %w", err)
			}
			if rec == nil {
				break
			}
			seqs = append(seqs, rec.FullSeq())
			names = append(names, rec.Name())
		}

		if len(seqs) == 0 {
			return fmt.Errorf("no sequences found in input file")
		}

		// Configure alignment options
		alignOpts := align.OntAlignmentDefaults()
		if !msaONT {
			alignOpts = align.DnaAlignmentDefaults()
		}
		alignOpts = alignOpts.Verbose(msaVerbose)

		opts := align.NewMSAOptions(alignOpts).
			MaxWorkers(msaThreads).
			Verbose(msaVerbose)

		profile := align.MSA(seqs, opts)
		if profile == nil {
			return fmt.Errorf("MSA produced no result")
		}

		// Open output
		var out *os.File
		if msaOutput == "" || msaOutput == "-" {
			out = os.Stdout
		} else {
			out, err = os.Create(msaOutput)
			if err != nil {
				return fmt.Errorf("opening output: %w", err)
			}
			defer out.Close()
		}

		if msaConsensus {
			// Output consensus as a single FASTA record
			cons := profile.Consensus()
			fmt.Fprintf(out, ">consensus\n%s\n", cons)
		} else {
			// Output gapped MSA as multi-sequence FASTA
			gapped := profile.GappedSequences()
			for i, seq := range gapped {
				fmt.Fprintf(out, ">%s\n%s\n", profile.Names[i], seq)
			}
		}

		return nil
	},
}

var (
	msaONT       bool
	msaThreads   int
	msaOutput    string
	msaConsensus bool
	msaVerbose   bool
)

func init() {
	msaCmd.Flags().BoolVar(&msaONT, "ont", true, "Use Oxford Nanopore alignment defaults (set --ont=false for Illumina)")
	msaCmd.Flags().IntVarP(&msaThreads, "threads", "t", 1, "Max parallel workers for all-pairs alignment")
	msaCmd.Flags().StringVarP(&msaOutput, "output", "o", "", "Output file (default: stdout)")
	msaCmd.Flags().BoolVar(&msaConsensus, "consensus", false, "Output a single consensus sequence instead of the full MSA")
	msaCmd.Flags().BoolVarP(&msaVerbose, "verbose", "v", false, "Enable verbose output")
}
