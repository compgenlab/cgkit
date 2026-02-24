package cmd

import (
	_ "embed"
	"fmt"
	"strings"

	"github.com/compgen-io/cgltk/align"
	"github.com/compgen-io/cgltk/seqio"
	"github.com/compgen-io/cgltk/sequtils"
	"github.com/spf13/cobra"
)

// fastagcCmd implements the initial counting entrypoint.
var ontPrimersCmd = &cobra.Command{
	Use:   "ont-primers <input.fastq>",
	Short: "Find and trim common ONT primers from the start of reads in a FASTQ file",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		umiSeq := "TTVVVVTTVVVVTTVVVVTTVVVVTTTGGG"
		var primers = make(map[string]string)

		var ontPrimers seqio.SeqReader
		var err error
		if ontPrimersFilename != "" {
			ontPrimers, err = seqio.NewFastaFile(ontPrimersFilename)
		} else {
			ontPrimers, err = seqio.NewFastaReader(strings.NewReader(ontPrimersDefault))
		}
		if err != nil {
			return err
		}
		sspName := "SSP"
		if ontPrimersUMI {
			sspName = "SSPII"
		}
		for seq, err := ontPrimers.NextSeq(); err == nil && seq != nil; seq, err = ontPrimers.NextSeq() {
			primers[seq.Name()] = seq.FullSeq().Seq()
		}

		fqReader, err := seqio.NewFastqFile(args[0])

		vnpseq := primers["VNP"]
		vnprc := sequtils.ReverseCompliment(vnpseq)
		sspseq := primers[sspName]
		ssprc := sequtils.ReverseCompliment(sspseq)

		opts := align.OntAlignmentDefaults().ClippingDisable()
		aligner := align.NewLocalAligner(opts)

		fmt.Printf("read\tlength\tVNP_score\tVNP_start_end\t%s\tVNP_strand\tSSP_score\tSSP_start_end\t%s\tSSP_strand", vnpseq, sspseq)
		if ontPrimersUMI {
			fmt.Print("\tUMI\tUMI_code")
		}
		if ontPrimersBC {
			fmt.Print("\tbarcode\tbarcode_seq\tbarcode_score")
		}
		fmt.Println()

		for seq, err := fqReader.NextSeq(); err == nil && seq != nil; seq, err = fqReader.NextSeq() {
			vnpAln := aligner.Align(vnpseq, seq.FullSeq().Seq(), "VNP", seq.Name())
			vnpAlnRc := aligner.Align(vnprc, seq.FullSeq().Seq(), "VNP", seq.Name())
			vnpAlnRc.QueryRevComp = true

			sspAln := aligner.Align(sspseq, seq.FullSeq().Seq(), sspName, seq.Name())
			sspAlnRc := aligner.Align(ssprc, seq.FullSeq().Seq(), sspName, seq.Name())
			sspAlnRc.QueryRevComp = true

			vtSeq := vnpAln.TargetStr()
			vtStrand := "+"
			if vnpAlnRc.Score > vnpAln.Score {
				vnpAln = vnpAlnRc
				vtSeq = sequtils.ReverseCompliment(vnpAlnRc.TargetStr())
				vtStrand = "-"
			}

			stSeq := sspAln.TargetStr()
			stStrand := "+"
			if sspAlnRc.Score > sspAln.Score {
				sspAln = sspAlnRc
				stSeq = sequtils.ReverseCompliment(sspAlnRc.TargetStr())
				stStrand = "-"
			}

			fmt.Printf("@%s\t%d\t%.2g\t%d-%d\t%s\t%s\t%.2g\t%d-%d\t%s\t%s", seq.Name(), len(seq.FullSeq().Seq()), vnpAln.Score, vnpAln.TargetStart, vnpAln.TargetEnd, vtSeq, vtStrand,
				sspAln.Score, sspAln.TargetStart, sspAln.TargetEnd, stSeq, stStrand)

			if ontPrimersUMI {
				umiCode := ""
				umiAln := aligner.Align(umiSeq, stSeq, "UMI", "target-SSP")
				var umiStr strings.Builder
				qPos := umiAln.QueryStart
				tPos := umiAln.TargetStart
				cigarExpanded, err := align.CigarExpand(umiAln.CIGAR)
				if err == nil {
					for i := 0; i < len(cigarExpanded); i++ {
						// fmt.Printf("qStr: %s\ntStr: %s\n-\n", qStr, tStr)
						op := cigarExpanded[i]
						switch op {
						case 'M':
							umiStr.WriteString(string(umiAln.TargetSeq[tPos]))
							if umiAln.QuerySeq[qPos] == 'V' {
								umiCode += string(umiAln.TargetSeq[tPos])
							} else if umiCode != "" && umiCode[len(umiCode)-1] != ':' {
								umiCode += ":"
							}
							qPos++
							tPos++
						case 'D':
							umiStr.WriteString(string(umiAln.TargetSeq[tPos]))
							tPos++
						case 'I':
							qPos++
						}
					}
				}
				if umiCode[len(umiCode)-1] == ':' {
					umiCode = umiCode[0 : len(umiCode)-1]
				}
				fmt.Printf("\t%s\t%s", umiStr.String(), umiCode)
			}

			if ontPrimersBC {
				var flankseq string
				if !vnpAln.QueryRevComp {
					start := max(0, vnpAln.TargetStart-30)
					flankseq = vnpAln.TargetSeq[start:vnpAln.TargetStart]
				} else {
					end := min(len(vnpAln.TargetSeq), vnpAln.TargetEnd+30)
					flankseq = vnpAln.TargetSeq[vnpAln.TargetEnd:end]
					flankseq = sequtils.ReverseCompliment(flankseq)
				}
				var bestAln *align.PairwiseAlignment
				var bestBC string
				for k, v := range primers {
					if k == "VNP" || k == "SSP" || k == "SSPII" {
						continue
					}

					aln := aligner.Align(flankseq, v, "flanking", k)
					if bestAln == nil || aln.Score > bestAln.Score {
						bestAln = aln
						bestBC = k
					}
				}
				fmt.Printf("\t%s\t%s\t%.2g", bestBC, bestAln.QuerySeq[bestAln.QueryStart:bestAln.QueryEnd], bestAln.Score)
			}

			fmt.Println()
		}

		return nil
	},
}

var ontPrimersFilename string
var ontPrimersUMI bool
var ontPrimersBC bool

//go:embed data/ont_seq.fa
var ontPrimersDefault string

func init() {
	ontPrimersCmd.Flags().StringVar(&ontPrimersFilename, "fasta", "", "FASTA file with primers (default use included primers)")
	ontPrimersCmd.Flags().BoolVar(&ontPrimersUMI, "umi", false, "Use UMI SSP primer (SSPII)")
	ontPrimersCmd.Flags().BoolVar(&ontPrimersBC, "barcode", false, "Identify the barcode used (upstream of VNP)")
	rootCmd.AddCommand(ontPrimersCmd)
}
