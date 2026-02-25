package ontcmd

import (
	_ "embed"
	"fmt"
	"runtime"
	"strings"

	"github.com/compgen-io/cgltk/align"
	"github.com/compgen-io/cgltk/seqio"
	"github.com/compgen-io/cgltk/support/utils"
	"github.com/spf13/cobra"
)

// fastagcCmd implements the initial counting entrypoint.
var ontPrimersCmd = &cobra.Command{
	GroupID: "ontcmd",
	Use:     "ont-primers <input.fastq>",
	Short:   "Find and trim common ONT primers from the start of reads in a FASTQ file",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}

		if ontThreads == 0 {
			ontThreads = runtime.GOMAXPROCS(0)
		}

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
		var vnpseq, sspseq seqio.SeqQual
		var barcodeSeqs = make([]seqio.SeqQual, 0)
		umiSeq := seqio.NewStringSeqQual("TTVVVVTTVVVVTTVVVVTTVVVVTTTGGG", "", "UMI").FullSeq()

		for seq, err := ontPrimers.NextSeq(); err == nil && seq != nil; seq, err = ontPrimers.NextSeq() {
			switch seq.Name() {
			case "VNP":
				vnpseq = seq.FullSeq()
			case "SSP":
				if !ontPrimersUMI {
					sspseq = seq.FullSeq()
				}
			case "SSPII":
				if ontPrimersUMI {
					sspseq = seq.FullSeq()
				}
			default:
				barcodeSeqs = append(barcodeSeqs, seq.FullSeq())
			}
		}

		fqReader, err := seqio.NewFastqFile(args[0])

		opts := align.OntAlignmentDefaults().ClippingDisable()
		aligner := align.NewLocalAligner(opts)

		fmt.Printf("read\tlength\tVNP_score\tVNP_start_end\t%s\tVNP_strand\tSSP_score\tSSP_start_end\t%s\tSSP_strand", vnpseq.Seq(), sspseq.Seq())
		if ontPrimersUMI {
			fmt.Print("\tUMI\tUMI_code\tUMI_score")
		}
		if ontPrimersBC {
			fmt.Print("\tbarcode\tbarcode_seq\tbarcode_score")
		}
		fmt.Println()

		// this will act as a processor limit for batch processing
		sem := utils.NewSemaphore(ontThreads)

		for seq, err := fqReader.NextSeq(); err == nil && seq != nil; seq, err = fqReader.NextSeq() {
			vnpAlnPromise := align.AlignBatch(aligner, sem, []seqio.SeqQual{vnpseq}, []seqio.SeqQual{seq.FullSeq(), seq.FullSeq().RevComp()})
			sspAlnPromise := align.AlignBatch(aligner, sem, []seqio.SeqQual{sspseq}, []seqio.SeqQual{seq.FullSeq(), seq.FullSeq().RevComp()})

			vnpAln := vnpAlnPromise.Result()
			sspAln := sspAlnPromise.Result()

			fmt.Printf("@%s\t%d\t%.2g\t%d-%d\t%s\t%s\t%.2g\t%d-%d\t%s\t%s", seq.Name(), len(seq.FullSeq().Seq()),
				vnpAln.Score, vnpAln.TargetStart, vnpAln.TargetEnd,
				vnpAln.TargetStr(), vnpAln.Target.Strand(),
				sspAln.Score, sspAln.TargetStart, sspAln.TargetEnd, sspAln.TargetStr(), sspAln.Target.Strand())

			if ontPrimersUMI {
				umiCode := ""
				umiAln := aligner.Align(umiSeq, sspAln.TargetSub())
				qStr := umiAln.QueryAlignedStr()
				tStr := umiAln.TargetAlignedStr()

				vMin := 0
				vMax := len(qStr)

				for i := 0; i < len(qStr) && qStr[i] != 'V'; i++ {
					vMin = i
				}
				for i := len(qStr) - 1; i >= 0 && qStr[i] != 'V'; i-- {
					vMax = i
				}

				inT := false
				for i := vMin + 1; i < vMax; i++ {
					switch tStr[i] {
					case 'T':
						if !inT {
							umiCode += ":"
						}
						inT = true
					case '-':
					default:
						inT = false
						umiCode += string(tStr[i])
					}
				}
				fmt.Printf("\t%s\t%s\t%.2g", umiAln.TargetStr(), umiCode, umiAln.Score)
			}

			if ontPrimersBC {
				// the vnpAln.Target will always be relative to the VNP seq, so we should look 30bp upstream always
				// if this is on the revcomp strand, we still look 30bp up.
				start := max(0, vnpAln.TargetStart-30)
				flankseq := vnpAln.Target.Sub(start, vnpAln.TargetStart)

				bestBC := align.AlignBatch(aligner, sem, barcodeSeqs, []seqio.SeqQual{flankseq}).Result()
				if bestBC.Score > 0 {
					fmt.Printf("\t%s\t%s\t%.2g", bestBC.Query.Name(), bestBC.TargetStr(), bestBC.Score)
				} else {
					fmt.Print("\t\t\t")
				}
			}

			fmt.Println()
		}

		return nil
	},
}

var ontPrimersFilename string
var ontPrimersUMI bool
var ontPrimersBC bool

var ontThreads int

//go:embed data/ont_seq.fa
var ontPrimersDefault string

func init() {
	ontPrimersCmd.Flags().StringVar(&ontPrimersFilename, "fasta", "", "FASTA file with primers (default use included primers)")
	ontPrimersCmd.Flags().BoolVar(&ontPrimersUMI, "umi", false, "Use UMI SSP primer (SSPII)")
	ontPrimersCmd.Flags().BoolVar(&ontPrimersBC, "barcode", false, "Identify the barcode used (upstream of VNP)")
	ontPrimersCmd.Flags().IntVarP(&ontThreads, "threads", "t", 0, "Threads to use (default: CPU count)")
}
