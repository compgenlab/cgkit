package seqcmd

import (
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/compgenlab/cghts/align"
	"github.com/compgenlab/cghts/seqio"
	"github.com/spf13/cobra"
)

// seq-msa is a thin CLI wrapper around align.MSA. All of the MSA logic —
// HP compression, reference handling, CLUSTAL/FASTA output, rehydration —
// lives in the align package as library methods so that future commands
// can reuse the exact same pipeline.

var msaCmd = &cobra.Command{
	GroupID:     "seqcmd",
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "seq-msa <input.fasta|fastq>",
	Short:       "Multiple sequence alignment via incremental consensus",
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

		// --hp-expand is a post-processing step on top of --hp-compress.
		// Expansion needs the per-sequence HP length data that compression
		// records, so it cannot be used on its own.
		if msaHPExpand && !msaHPCompress {
			return fmt.Errorf("--hp-expand requires --hp-compress")
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

		// --hp-expand post-processes the alignment by expanding each
		// compressed column out to the per-row original HP lengths. When
		// --consensus is also set, the consensus is computed first and
		// included as an extra row in the expanded alignment (so the
		// per-column max covers the consensus runs too). After expansion
		// the alignment behaves like a normal full-length one, so the
		// rest of the output dispatch below doesn't have to branch on
		// HP mode separately.
		//
		// A subtle point: with --hp-expand --consensus, the consensus
		// row is embedded in the expanded alignment. So for FASTA output
		// (--fasta), the consensus row would be written as one of the
		// sequence records. To keep --fasta --consensus meaning "just
		// the consensus as a single record", we detect that combination
		// separately below and produce the single-record output without
		// invoking Expanded(withConsensus=true).
		displayAln := aln
		if msaHPExpand && !(msaFasta && msaConsensus) {
			// Expand the alignment in place. Include the consensus row
			// only when we're about to show it (CLUSTAL + --consensus).
			displayAln = aln.Expanded(msaConsensus)
			if displayAln == nil {
				return fmt.Errorf("Expanded() returned nil; did --hp-compress populate HPLens?")
			}
		}

		// Output dispatch — --consensus is a modifier on top of --fasta
		// or the default CLUSTAL output, not a standalone mode:
		//
		//   default             : CLUSTAL interleaved (compressed or
		//                          expanded depending on --hp-expand).
		//   --consensus         : CLUSTAL + a synthetic consensus row
		//                          at the bottom of every block. With
		//                          --hp-expand the consensus is already
		//                          embedded in displayAln as a regular
		//                          row, so we fall through to plain
		//                          WriteClustal.
		//   --fasta             : gapped multi-FASTA of the alignment
		//                          rows.
		//   --fasta --consensus : single-record FASTA of the consensus
		//                          sequence. Rehydrated when --hp-expand
		//                          (or just --hp-compress) is set.
		switch {
		case msaFasta && msaConsensus:
			var cons string
			if msaHPCompress {
				cons = aln.RehydratedConsensus()
			} else {
				cons = aln.Consensus()
			}
			fmt.Fprintf(out, ">consensus\n%s\n", cons)
		case msaFasta:
			if err := displayAln.WriteFasta(out); err != nil {
				return fmt.Errorf("writing fasta: %w", err)
			}
		case msaConsensus:
			if msaHPExpand {
				// Consensus is already a row in displayAln after
				// Expanded(true); emit it with plain WriteClustal so we
				// don't get a duplicate synthetic consensus row.
				if err := displayAln.WriteClustal(out); err != nil {
					return fmt.Errorf("writing clustal: %w", err)
				}
			} else {
				if err := displayAln.WriteClustalWithConsensus(out); err != nil {
					return fmt.Errorf("writing clustal: %w", err)
				}
			}
		default:
			if err := displayAln.WriteClustal(out); err != nil {
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
	msaHPExpand   bool
	msaRef        string
	msaVerbose    bool
)

func init() {
	msaCmd.Flags().BoolVar(&msaONT, "ont", true, "Use Oxford Nanopore alignment defaults (set --ont=false for Illumina)")
	msaCmd.Flags().IntVarP(&msaThreads, "threads", "t", 1, "Max parallel workers for all-pairs alignment")
	msaCmd.Flags().StringVarP(&msaOutput, "output", "o", "", "Output file (default: stdout)")
	msaCmd.Flags().BoolVar(&msaFasta, "fasta", false, "Output gapped multi-sequence FASTA instead of CLUSTAL")
	msaCmd.Flags().BoolVar(&msaConsensus, "consensus", false, "Append a consensus row to the CLUSTAL output (or with --fasta, write only the consensus as FASTA)")
	msaCmd.Flags().BoolVar(&msaHPCompress, "hp-compress", false, "Homopolymer-compress sequences before alignment")
	msaCmd.Flags().BoolVar(&msaHPExpand, "hp-expand", false, "After HP-compressed MSA, expand each row back to its original homopolymer lengths (requires --hp-compress)")
	msaCmd.Flags().StringVar(&msaRef, "ref", "", "Name of the reference sequence in the input; aligned last and shown first, used for HP tiebreaks")
	msaCmd.Flags().BoolVarP(&msaVerbose, "verbose", "v", false, "Enable verbose output")
}
