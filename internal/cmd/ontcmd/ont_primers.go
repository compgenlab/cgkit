package ontcmd

import (
	_ "embed"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"runtime"
	"strings"

	"github.com/compgen-io/cgltk/align"
	"github.com/compgen-io/cgltk/seqio"
	"github.com/compgen-io/cgltk/support/utils"
	"github.com/spf13/cobra"
)

// openWriter opens a file for writing (or returns os.Stdout if filename is empty).
// If the filename ends in ".gz", the output is gzip-compressed.
// The returned closer must be called when done.
func openWriter(filename string) (io.Writer, func() error, error) {
	if filename == "" {
		return nil, func() error { return nil }, nil
	}
	if filename == "-" {
		return os.Stdout, func() error { return nil }, nil
	}
	f, err := os.Create(filename)
	if err != nil {
		return nil, nil, err
	}
	if strings.HasSuffix(filename, ".gz") {
		gz := gzip.NewWriter(f)
		return gz, func() error {
			if err := gz.Close(); err != nil {
				f.Close()
				return err
			}
			return f.Close()
		}, nil
	}
	return f, f.Close, nil
}

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
		// --add-umi implies --detect-umi; --add-barcode implies --detect-barcode
		if ontWriteUMI {
			ontPrimersUMI = true
		}
		if ontWriteBarcode {
			ontDetectBC = true
		}

		var vnpseq, sspseq seqio.SeqQual
		var barcodeSeqs = make([]seqio.SeqQual, 0)
		umiSeq := seqio.NewStringSeq("TTVVVVTTVVVVTTVVVVTTVVVVTTTGGG", "UMI").FullSeq()

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

		if vnpseq.Seq() == "" {
			return fmt.Errorf("VNP sequence not found in primer FASTA")
		}
		if sspseq.Seq() == "" {
			if ontPrimersUMI {
				return fmt.Errorf("SSPII sequence not found in primer FASTA")
			}
			return fmt.Errorf("SSP sequence not found in primer FASTA")
		}

		// Build accepted barcode set from comma-list; non-empty list implies detection.
		if len(barcodeSeqs) == 0 && (ontDetectBC || ontFilterBarcodes != "" || ontFilterBarcodeScore >= 0 || ontFilterBarcodeMatches >= 0 || ontWriteBarcode) {
			return fmt.Errorf("barcode detection requested but no barcode sequences found in primer FASTA")
		}

		acceptedBarcodes := make(map[string]bool)
		for _, bc := range strings.Split(ontFilterBarcodes, ",") {
			bc = strings.TrimSpace(bc)
			if bc != "" {
				acceptedBarcodes[bc] = true
			}
		}
		needsBarcode := ontDetectBC || len(acceptedBarcodes) > 0 ||
			ontFilterBarcodeScore >= 0 || ontFilterBarcodeMatches >= 0 || ontWriteBarcode

		// Open report writer (default: stdout).
		reportWriter, closeReport, err := openWriter(ontReportFilename)
		if err != nil {
			return fmt.Errorf("opening report: %w", err)
		}
		defer func() {
			if err := closeReport(); err != nil {
				fmt.Fprintf(os.Stderr, "error closing report: %v\n", err)
			}
		}()

		// Open passing/failed FASTQ writers.
		var passingWriter *seqio.FastqWriter
		if ontPassingFQFilename != "" {
			passingWriter, err = seqio.OpenFastqWriter(ontPassingFQFilename)
			if err != nil {
				return fmt.Errorf("opening passing-fastq: %w", err)
			}
			defer func() {
				if err := passingWriter.Close(); err != nil {
					fmt.Fprintf(os.Stderr, "error closing passing-fastq: %v\n", err)
				}
			}()
		}
		var failedWriter *seqio.FastqWriter
		if ontFailedFQFilename != "" {
			failedWriter, err = seqio.OpenFastqWriter(ontFailedFQFilename)
			if err != nil {
				return fmt.Errorf("opening failed-fastq: %w", err)
			}
			defer func() {
				if err := failedWriter.Close(); err != nil {
					fmt.Fprintf(os.Stderr, "error closing failed-fastq: %v\n", err)
				}
			}()
		}

		if reportWriter == nil && passingWriter == nil && failedWriter == nil {
			return fmt.Errorf("no output specified: at least one of --report, --passing-fastq, or --failed-fastq is required")
		}

		if passingWriter == nil && failedWriter == nil {
			if ontWriteUMI {
				fmt.Fprintln(os.Stderr, "warning: --add-umi has no effect without --passing-fastq or --failed-fastq")
			}
			if ontWriteBarcode {
				fmt.Fprintln(os.Stderr, "warning: --add-barcode has no effect without --passing-fastq or --failed-fastq")
			}
			if ontTrimFlanking {
				fmt.Fprintln(os.Stderr, "warning: --trim-flanking has no effect without --passing-fastq")
			}
		} else if passingWriter == nil && ontTrimFlanking {
			fmt.Fprintln(os.Stderr, "warning: --trim-flanking has no effect without --passing-fastq")
		}

		fqReader, err := seqio.NewFastqFile(args[0])
		if err != nil {
			return err
		}

		opts := align.OntAlignmentDefaults().ClippingDisable()
		aligner := align.NewLocalAligner(opts)

		if reportWriter != nil {
			fmt.Fprintf(reportWriter, "read\tlength\tVNP_score\tVNP_matches\tVNP_start\tVNP_end\t%s\tVNP_strand\tSSP_score\tSSP_matches\tSSP_start\tSSP_end\t%s\tSSP_strand", vnpseq.Seq(), sspseq.Seq())
			if ontPrimersUMI {
				fmt.Fprint(reportWriter, "\tUMI\tUMI_code\tUMI_score")
			}
			if needsBarcode {
				fmt.Fprint(reportWriter, "\tbarcode\tbarcode_seq\tbarcode_score\tbarcode_matches")
			}
			fmt.Fprintln(reportWriter)
		}

		sem := utils.NewSemaphore(ontThreads)

		for seq, err := fqReader.NextFastqSeq(); err == nil && seq != nil; seq, err = fqReader.NextFastqSeq() {
			seqFull := seq.FullSeq()
			seqLen := len(seqFull.Seq())
			seqRevComp := seqFull.RevComp()
			vnpAlnPromise := align.AlignBatch(aligner, sem, []seqio.SeqQual{vnpseq}, []seqio.SeqQual{seqFull, seqRevComp})
			sspAlnPromise := align.AlignBatch(aligner, sem, []seqio.SeqQual{sspseq}, []seqio.SeqQual{seqFull, seqRevComp})

			vnpAln := vnpAlnPromise.Result()
			sspAln := sspAlnPromise.Result()

			vnpStart, vnpEnd := vnpAln.TargetStart, vnpAln.TargetEnd
			if vnpAln.Target.IsRevComp() {
				vnpStart = seqLen - vnpStart
				vnpEnd = seqLen - vnpEnd
				vnpStart, vnpEnd = vnpEnd, vnpStart
			}
			sspStart, sspEnd := sspAln.TargetStart, sspAln.TargetEnd
			if sspAln.Target.IsRevComp() {
				sspStart = seqLen - sspStart
				sspEnd = seqLen - sspEnd
				sspStart, sspEnd = sspEnd, sspStart
			}

			if reportWriter != nil {
				fmt.Fprintf(reportWriter, "@%s\t%d\t%.2g\t%d\t%d\t%d\t%s\t%s\t%.2g\t%d\t%d\t%d\t%s\t%s",
					seq.Name(), seqLen,
					vnpAln.Score, vnpAln.Matches(), vnpStart, vnpEnd,
					vnpAln.TargetStr(), vnpAln.Target.Strand(),
					sspAln.Score, sspAln.Matches(), sspStart, sspEnd,
					sspAln.TargetStr(), sspAln.Target.Strand())
			}

			// UMI processing.
			var umiCode string
			if ontPrimersUMI {
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

				inT := true
				for i := vMin + 1; i < vMax; i++ {
					switch tStr[i] {
					case 'T':
						if !inT {
							umiCode += "-"
						}
						inT = true
					case '-':
					default:
						inT = false
						umiCode += string(tStr[i])
					}
				}
				if reportWriter != nil {
					fmt.Fprintf(reportWriter, "\t%s\t%s\t%.2g", umiAln.TargetStr(), umiCode, umiAln.Score)
				}
			}

			// Barcode processing.
			var bestBCName, bestBCSeq string
			var bestBCScore float32
			var bestBCMatches int
			if needsBarcode {
				start := max(0, vnpAln.TargetStart-32)
				flankseq := vnpAln.Target.Sub(start, vnpAln.TargetStart)
				bestBC := align.AlignBatch(aligner, sem, barcodeSeqs, []seqio.SeqQual{flankseq}).Result()
				if bestBC.Score > 0 {
					bestBCName = bestBC.Query.Name()
					bestBCSeq = bestBC.TargetStr()
					bestBCScore = bestBC.Score
					bestBCMatches = bestBC.Matches()
					if reportWriter != nil {
						fmt.Fprintf(reportWriter, "\t%s\t%s\t%.2g\t%d", bestBCName, bestBCSeq, bestBCScore, bestBCMatches)
					}
				} else if reportWriter != nil {
					fmt.Fprint(reportWriter, "\t\t\t\t")
				}
			}
			if reportWriter != nil {
				fmt.Fprintln(reportWriter)
			}

			// Write to passing or failed FASTQ if requested.
			if passingWriter != nil || failedWriter != nil {
				passing := true
				var failReasons []string

				if ontFilterVNPMatches >= 0 && vnpAln.Matches() < ontFilterVNPMatches {
					passing = false
					failReasons = append(failReasons, "low_vnp_match")
				}
				if ontFilterSSPMatches >= 0 && sspAln.Matches() < ontFilterSSPMatches {
					passing = false
					failReasons = append(failReasons, "low_ssp_match")
				}
				if ontFilterVNPScore >= 0 && vnpAln.Score < ontFilterVNPScore {
					passing = false
					failReasons = append(failReasons, "low_vnp_score")
				}
				if ontFilterSSPScore >= 0 && sspAln.Score < ontFilterSSPScore {
					passing = false
					failReasons = append(failReasons, "low_ssp_score")
				}
				if ontFilterVNPSSPPair && vnpAln.Target.Strand() == sspAln.Target.Strand() {
					passing = false // valid pairs flank on opposite strands
					failReasons = append(failReasons, "unpaired")
				}
				if len(acceptedBarcodes) > 0 && !acceptedBarcodes[bestBCName] {
					passing = false
					failReasons = append(failReasons, "barcode_mismatch")
				}
				if ontFilterBarcodeScore >= 0 && bestBCScore < ontFilterBarcodeScore {
					passing = false
					failReasons = append(failReasons, "low_barcode_score")
				}
				if ontFilterBarcodeMatches >= 0 && bestBCMatches < ontFilterBarcodeMatches {
					passing = false
					failReasons = append(failReasons, "low_barcode_match")
				}

				// Annotate FASTQ comment with SAM-style tags.
				if passing {
					if ontWriteUMI && umiCode != "" {
						seq.AddCommentTSV("RX:Z:" + umiCode)
					}
					if ontWriteBarcode && bestBCName != "" {
						seq.AddCommentTSV("BC:Z:" + bestBCSeq)
						seq.AddCommentTSV("ZB:Z:" + bestBCName)
					}
				} else {
					seq.AddCommentTSV("CO:Z:" + strings.Join(failReasons, ";"))
				}

				writer := failedWriter
				if passing {
					writer = passingWriter
				}
				if writer != nil {
					if passing && ontTrimFlanking {
						vnpValid := vnpAln.Score > 0 &&
							(ontFilterVNPScore < 0 || vnpAln.Score >= ontFilterVNPScore) &&
							(ontFilterVNPMatches < 0 || vnpAln.Matches() >= ontFilterVNPMatches)
						sspValid := sspAln.Score > 0 &&
							(ontFilterSSPScore < 0 || sspAln.Score >= ontFilterSSPScore) &&
							(ontFilterSSPMatches < 0 || sspAln.Matches() >= ontFilterSSPMatches)

						trimStart := 0
						trimEnd := seqLen

						vnpPlus := !vnpAln.Target.IsRevComp()
						sspPlus := !sspAln.Target.IsRevComp()

						if vnpValid && sspValid && vnpPlus != sspPlus {
							// Both found on opposite strands — trim both ends.
							if vnpPlus {
								trimStart = vnpEnd
								trimEnd = sspStart
							} else {
								trimStart = sspEnd
								trimEnd = vnpStart
							}
						} else if vnpValid && !sspValid {
							// Only VNP found.
							if vnpPlus {
								trimStart = vnpEnd
							} else {
								trimEnd = vnpStart
							}
						} else if sspValid && !vnpValid {
							// Only SSP found.
							if sspPlus {
								trimStart = sspEnd
							} else {
								trimEnd = sspStart
							}
						}
						// Both on same strand: no trimming (trimStart/trimEnd unchanged).

						if trimStart < trimEnd {
							trimmed := seqFull.Sub(trimStart, trimEnd)
							if err := writer.WriteRecord(seq.Name(), seq.Comment(), trimmed.Seq(), trimmed.Qual()); err != nil {
								return err
							}
						} else {
							// Trimming produced empty or inverted range — write as-is.
							if err := writer.Write(seq); err != nil {
								return err
							}
						}
					} else {
						if err := writer.Write(seq); err != nil {
							return err
						}
					}
				}
			}
		}

		return nil
	},
}

var ontPrimersFilename string
var ontPassingFQFilename string
var ontFailedFQFilename string
var ontReportFilename string
var ontPrimersUMI bool
var ontDetectBC bool
var ontFilterBarcodes string

var ontFilterVNPScore float32
var ontFilterSSPScore float32
var ontFilterBarcodeScore float32
var ontFilterVNPSSPPair bool

var ontFilterVNPMatches int
var ontFilterSSPMatches int
var ontFilterBarcodeMatches int

var ontWriteBarcode bool
var ontWriteUMI bool
var ontTrimFlanking bool

var ontThreads int

//go:embed data/ont_seq.fa
var ontPrimersDefault string

func init() {
	ontPrimersCmd.Flags().StringVar(&ontPassingFQFilename, "passing-fastq", "", "Write passing reads to this file (gzipped if .gz)")
	ontPrimersCmd.Flags().StringVar(&ontFailedFQFilename, "failed-fastq", "", "Write failed reads to this file (gzipped if .gz)")
	ontPrimersCmd.Flags().StringVar(&ontReportFilename, "report", "", "Write tab-delimited report to this file (use '-' for stdout; gzipped if .gz)")
	ontPrimersCmd.Flags().StringVar(&ontPrimersFilename, "primers-fasta", "", "FASTA file with primers (default: use included primers)")

	ontPrimersCmd.Flags().BoolVar(&ontWriteBarcode, "add-barcode", false, "Add BC= tag to FASTQ comment when writing output")
	ontPrimersCmd.Flags().BoolVar(&ontWriteUMI, "add-umi", false, "Add UMI= tag to FASTQ comment when writing output")
	ontPrimersCmd.Flags().BoolVar(&ontTrimFlanking, "trim-flanking", false, "Trim VNP/SSP sequences from passing reads before writing")

	ontPrimersCmd.Flags().BoolVar(&ontFilterVNPSSPPair, "filter-pair", false, "Require paired VNP/SSP (flanking on opposite strands)")
	ontPrimersCmd.Flags().Float32Var(&ontFilterVNPScore, "filter-vnp-score", -1, "Require minimum VNP alignment score")
	ontPrimersCmd.Flags().Float32Var(&ontFilterSSPScore, "filter-ssp-score", -1, "Require minimum SSP alignment score")
	ontPrimersCmd.Flags().Float32Var(&ontFilterBarcodeScore, "filter-barcode-score", -1, "Require minimum barcode alignment score")

	ontPrimersCmd.Flags().IntVar(&ontFilterVNPMatches, "filter-vnp-match", -1, "Require minimum VNP match count")
	ontPrimersCmd.Flags().IntVar(&ontFilterSSPMatches, "filter-ssp-match", -1, "Require minimum SSP match count")
	ontPrimersCmd.Flags().IntVar(&ontFilterBarcodeMatches, "filter-barcode-match", -1, "Require minimum barcode match count")

	ontPrimersCmd.Flags().StringVar(&ontFilterBarcodes, "filter-allowed-barcodes", "", "Comma-separated list of acceptable barcode names (also enables barcode detection)")
	ontPrimersCmd.Flags().BoolVar(&ontDetectBC, "detect-barcode", false, "Detect barcode upstream of VNP and include in report (no name filtering)")

	ontPrimersCmd.Flags().BoolVar(&ontPrimersUMI, "detect-umi", false, "Use UMI SSP primer (SSPII)")
	ontPrimersCmd.Flags().IntVarP(&ontThreads, "threads", "t", 0, "Threads to use (default: CPU count)")
}
