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

// seq-msa is a thin CLI wrapper around align.MSA. All of the MSA logic —
// HP compression, reference handling, CLUSTAL/FASTA output, rehydration —
// lives in the align package as library methods so that future commands
// can reuse the exact same pipeline.

var msaCmd = &cobra.Command{
	GroupID: "seqcmd",
	Use:     "seq-msa <input.fasta|fastq>",
	Short:   "Multiple sequence alignment via incremental consensus",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}

		// Detect file format from extension. FASTQ support uses the same
		// reader interface, so downstream code does not need to care.
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

		// Load all sequences into memory. The MSA algorithm is in-memory
		// anyway (all-pairs alignment), so there's no streaming benefit.
		var seqs []seqio.SeqQual
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
		}

		if len(seqs) == 0 {
			return fmt.Errorf("no sequences found in input file")
		}

		// Build alignment options: scoring preset, threading, HP compression,
		// and reference name are all configured here and passed to align.MSA.
		alignOpts := align.OntAlignmentDefaults()
		if !msaONT {
			alignOpts = align.DnaAlignmentDefaults()
		}
		alignOpts = alignOpts.Verbose(msaVerbose)

		opts := align.NewMSAOptions(alignOpts).
			MaxWorkers(msaThreads).
			Verbose(msaVerbose).
			HPCompress(msaHPCompress).
			RefName(msaRef)

		aln, err := align.MSA(seqs, opts)
		if err != nil {
			return err
		}
		if aln == nil {
			return fmt.Errorf("MSA produced no result")
		}

		// Open output — file or stdout.
		var out io.Writer = os.Stdout
		if msaOutput != "" && msaOutput != "-" {
			f, err := os.Create(msaOutput)
			if err != nil {
				return fmt.Errorf("opening output: %w", err)
			}
			defer f.Close()
			out = f
		}

		// Dispatch to the library's output methods. All three formats live
		// on MSAAlignment so any caller can reuse them.
		//
		//   --consensus : single FASTA record. In HP-compressed mode this
		//                 is the rehydrated consensus (HP runs restored
		//                 from per-column mode length); otherwise it is
		//                 the plain majority-vote consensus (reads only,
		//                 ref excluded when present).
		//   --fasta     : gapped multi-FASTA of the alignment rows as-is.
		//                 In HP-compressed mode this shows the compressed
		//                 alignment; rehydration only applies to consensus.
		//   default     : CLUSTAL interleaved format, same rules as --fasta
		//                 for HP-compressed output.
		switch {
		case msaConsensus:
			var cons string
			if msaHPCompress {
				cons = aln.RehydratedConsensus()
			} else {
				cons = aln.Consensus()
			}
			fmt.Fprintf(out, ">consensus\n%s\n", cons)
		case msaFasta:
			if err := aln.WriteFasta(out); err != nil {
				return fmt.Errorf("writing fasta: %w", err)
			}
		default:
			if err := aln.WriteClustal(out); err != nil {
				return fmt.Errorf("writing clustal: %w", err)
			}
		}

		return nil
	},
}

var (
	msaONT        bool
	msaThreads    int
	msaOutput     string
	msaFasta      bool
	msaConsensus  bool
	msaHPCompress bool
	msaRef        string
	msaVerbose    bool
)

func init() {
	msaCmd.Flags().BoolVar(&msaONT, "ont", true, "Use Oxford Nanopore alignment defaults (set --ont=false for Illumina)")
	msaCmd.Flags().IntVarP(&msaThreads, "threads", "t", 1, "Max parallel workers for all-pairs alignment")
	msaCmd.Flags().StringVarP(&msaOutput, "output", "o", "", "Output file (default: stdout)")
	msaCmd.Flags().BoolVar(&msaFasta, "fasta", false, "Output gapped multi-sequence FASTA instead of CLUSTAL")
	msaCmd.Flags().BoolVar(&msaConsensus, "consensus", false, "Output a single consensus sequence instead of the full MSA")
	msaCmd.Flags().BoolVar(&msaHPCompress, "hp-compress", false, "Homopolymer-compress sequences before alignment; --consensus output is rehydrated using per-column mode lengths")
	msaCmd.Flags().StringVar(&msaRef, "ref", "", "Name of the reference sequence in the input; aligned last and shown first, used for HP tiebreaks")
	msaCmd.Flags().BoolVarP(&msaVerbose, "verbose", "v", false, "Enable verbose output")
}
