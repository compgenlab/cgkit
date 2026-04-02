package ontcmd

import (
	"compress/gzip"
	_ "embed"
	"fmt"
	"io"
	"os"
	"runtime"
	"strings"
	"sync"

	"github.com/compgen-io/cgltk/align"
	"github.com/compgen-io/cgltk/seqio"
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

		// Validate that every requested barcode name exists in the FASTA file.
		if len(acceptedBarcodes) > 0 {
			bcNames := make(map[string]bool, len(barcodeSeqs))
			for _, bc := range barcodeSeqs {
				bcNames[bc.Name()] = true
			}
			for name := range acceptedBarcodes {
				if !bcNames[name] {
					return fmt.Errorf("barcode %q not found in primer FASTA", name)
				}
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
			fmt.Fprint(reportWriter, "\tinsert_length\tstatus")
			fmt.Fprintln(reportWriter)
		}

		// Per-read result type, computed by workers.
		type readResult struct {
			seq           *seqio.FastqSeqRecord
			seqFull       seqio.SeqQual
			seqLen        int
			vnpAln        *align.PairwiseAlignment
			sspAln        *align.PairwiseAlignment
			vnpStart      int
			vnpEnd        int
			sspStart      int
			sspEnd        int
			insertLength  int
			umiTargetStr  string
			umiCode       string
			umiScore      float32
			bestBCName    string
			bestBCSeq     string
			bestBCScore   float32
			bestBCMatches int
		}

		type workItem struct {
			seq      *seqio.FastqSeqRecord
			resultCh chan *readResult
		}

		workCh := make(chan workItem, ontThreads)
		orderCh := make(chan chan *readResult, ontThreads)

		// Start worker goroutines. The aligner is thread-safe (all mutable state
		// is stack-local per Align() call), so workers share a single instance.
		var workerWg sync.WaitGroup
		for range ontThreads {
			workerWg.Add(1)
			go func() {
				defer workerWg.Done()
				for item := range workCh {
					seq := item.seq
					seqFull := seq.FullSeq()
					seqLen := len(seqFull.Seq())
					seqRevComp := seqFull.RevComp()

					// VNP: align against both strands, keep best.
					vnpAlnPlus := aligner.Align(vnpseq, seqFull)
					vnpAlnMinus := aligner.Align(vnpseq, seqRevComp)
					vnpAln := vnpAlnPlus
					if vnpAlnMinus.Score > vnpAlnPlus.Score {
						vnpAln = vnpAlnMinus
					}

					// SSP: align against both strands, keep best.
					sspAlnPlus := aligner.Align(sspseq, seqFull)
					sspAlnMinus := aligner.Align(sspseq, seqRevComp)
					sspAln := sspAlnPlus
					if sspAlnMinus.Score > sspAlnPlus.Score {
						sspAln = sspAlnMinus
					}

					// Convert RevComp target positions to plus-strand coordinates.
					vnpStart, vnpEnd := vnpAln.TargetStart, vnpAln.TargetEnd
					if vnpAln.Target.IsRevComp() {
						vnpStart = seqLen - vnpAln.TargetEnd
						vnpEnd = seqLen - vnpAln.TargetStart
					}
					sspStart, sspEnd := sspAln.TargetStart, sspAln.TargetEnd
					if sspAln.Target.IsRevComp() {
						sspStart = seqLen - sspAln.TargetEnd
						sspEnd = seqLen - sspAln.TargetStart
					}

					// Compute insert length: bases between VNP end and SSP start
					// when both primers flank on opposite strands.
					insertLength := -1
					vnpPlus := !vnpAln.Target.IsRevComp()
					sspPlus := !sspAln.Target.IsRevComp()
					if vnpPlus != sspPlus {
						if vnpPlus {
							insertLength = sspStart - vnpEnd
						} else {
							insertLength = vnpStart - sspEnd
						}
						if insertLength < 0 {
							insertLength = 0
						}
					}

					result := &readResult{
						seq:          seq,
						seqFull:      seqFull,
						seqLen:       seqLen,
						vnpAln:       vnpAln,
						sspAln:       sspAln,
						vnpStart:     vnpStart,
						vnpEnd:       vnpEnd,
						sspStart:     sspStart,
						sspEnd:       sspEnd,
						insertLength: insertLength,
					}

					// UMI processing.
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
									result.umiCode += "-"
								}
								inT = true
							case '-':
							default:
								inT = false
								result.umiCode += string(tStr[i])
							}
						}
						if ontUMISepT {
							result.umiCode = strings.ReplaceAll(result.umiCode, "-", "TT")
						}
						result.umiTargetStr = umiAln.TargetStr()
						result.umiScore = umiAln.Score
					}

					// Barcode processing.
					if needsBarcode {
						start := max(0, vnpAln.TargetStart-32)
						flankseq := vnpAln.Target.Sub(start, vnpAln.TargetStart)
						var bestBC *align.PairwiseAlignment
						for _, bc := range barcodeSeqs {
							aln := aligner.Align(bc, flankseq)
							if bestBC == nil || aln.Score > bestBC.Score {
								bestBC = aln
							}
						}
						if bestBC != nil && bestBC.Score > 0 {
							result.bestBCName = bestBC.Query.Name()
							result.bestBCSeq = bestBC.TargetStr()
							result.bestBCScore = bestBC.Score
							result.bestBCMatches = bestBC.Matches()
						}
					}

					item.resultCh <- result
				}
			}()
		}

		// Writer goroutine: ranges over orderCh in submission order, blocking on
		// each resultCh so output order matches input order.
		writerErrCh := make(chan error, 1)
		go func() {
			var writeErr error
			for resultCh := range orderCh {
				result := <-resultCh
				if writeErr != nil {
					continue // drain remaining results but skip writing
				}

				seq := result.seq
				vnpAln := result.vnpAln
				sspAln := result.sspAln

				// Compute pass/fail status (shared by report and FASTQ output).
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
				if len(acceptedBarcodes) > 0 && !acceptedBarcodes[result.bestBCName] {
					passing = false
					failReasons = append(failReasons, "barcode_mismatch")
				}
				if ontFilterBarcodeScore >= 0 && result.bestBCScore < ontFilterBarcodeScore {
					passing = false
					failReasons = append(failReasons, "low_barcode_score")
				}
				if ontFilterBarcodeMatches >= 0 && result.bestBCMatches < ontFilterBarcodeMatches {
					passing = false
					failReasons = append(failReasons, "low_barcode_match")
				}

				statusStr := "PASS"
				if !passing {
					statusStr = strings.Join(failReasons, ";")
				}

				if reportWriter != nil {
					fmt.Fprintf(reportWriter, "@%s\t%d\t%.2g\t%d\t%d\t%d\t%s\t%s\t%.2g\t%d\t%d\t%d\t%s\t%s",
						seq.Name(), result.seqLen,
						vnpAln.Score, vnpAln.Matches(), result.vnpStart, result.vnpEnd,
						vnpAln.TargetStr(), vnpAln.Target.Strand(),
						sspAln.Score, sspAln.Matches(), result.sspStart, result.sspEnd,
						sspAln.TargetStr(), sspAln.Target.Strand())

					if ontPrimersUMI {
						fmt.Fprintf(reportWriter, "\t%s\t%s\t%.2g", result.umiTargetStr, result.umiCode, result.umiScore)
					}

					if needsBarcode {
						if result.bestBCName != "" {
							fmt.Fprintf(reportWriter, "\t%s\t%s\t%.2g\t%d", result.bestBCName, result.bestBCSeq, result.bestBCScore, result.bestBCMatches)
						} else {
							fmt.Fprint(reportWriter, "\t\t\t\t")
						}
					}
					fmt.Fprintf(reportWriter, "\t%d\t%s\n", result.insertLength, statusStr)
				}

				// Write to passing or failed FASTQ if requested.
				if passingWriter != nil || failedWriter != nil {
					// Annotate FASTQ comment with SAM-style tags.
					if ontWriteUMI && result.umiCode != "" {
						seq.AddCommentTSV("RX:Z:" + result.umiCode)
					}
					if ontWriteBarcode && result.bestBCName != "" {
						seq.AddCommentTSV("BC:Z:" + result.bestBCSeq)
						seq.AddCommentTSV("ZB:Z:" + result.bestBCName)
					}
					if !passing {
						seq.AddCommentTSV("CO:Z:" + statusStr)
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
							trimEnd := result.seqLen

							vnpPlus := !vnpAln.Target.IsRevComp()
							sspPlus := !sspAln.Target.IsRevComp()

							if vnpValid && sspValid && vnpPlus != sspPlus {
								// Both found on opposite strands — trim both ends.
								if vnpPlus {
									trimStart = result.vnpEnd
									trimEnd = result.sspStart
								} else {
									trimStart = result.sspEnd
									trimEnd = result.vnpStart
								}
							} else if vnpValid && !sspValid {
								// Only VNP found.
								if vnpPlus {
									trimStart = result.vnpEnd
								} else {
									trimEnd = result.vnpStart
								}
							} else if sspValid && !vnpValid {
								// Only SSP found.
								if sspPlus {
									trimStart = result.sspEnd
								} else {
									trimEnd = result.sspStart
								}
							}
							// Both on same strand: no trimming (trimStart/trimEnd unchanged).

							if trimStart < trimEnd {
								trimmed := result.seqFull.Sub(trimStart, trimEnd)
								if ontSenseCorrect && vnpPlus {
									trimmed = trimmed.RevComp()
								}
								writeErr = writer.WriteRecord(seq.Name(), seq.Comment(), trimmed.Seq(), trimmed.Qual())
							} else {
								// Trimming produced empty or inverted range — write as-is.
								fullSeq := seq.FullSeq()
								if ontSenseCorrect && vnpPlus {
									fullSeq = seq.FullSeq().RevComp()
								}
								writeErr = writer.WriteRecord(seq.Name(), seq.Comment(), fullSeq.Seq(), fullSeq.Qual())
							}
						} else {
							writeErr = writer.Write(seq)
						}
					}
				}
			}
			writerErrCh <- writeErr
		}()

		// Reader: send each read to orderCh (for ordering) then workCh (for processing).
		// orderCh is sent first so the writer always sees results in input order.
		var readErr error
		for {
			seq, err := fqReader.NextFastqSeq()
			if err != nil {
				if err != io.EOF {
					readErr = err
				}
				break
			}
			resultCh := make(chan *readResult, 1)
			orderCh <- resultCh
			workCh <- workItem{seq: seq, resultCh: resultCh}
		}
		close(workCh)
		workerWg.Wait()
		close(orderCh)

		writeErr := <-writerErrCh
		if readErr != nil {
			return readErr
		}
		return writeErr
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
var ontUMISepT bool
var ontTrimFlanking bool
var ontSenseCorrect bool

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
	ontPrimersCmd.Flags().BoolVar(&ontUMISepT, "umi-sep-t", false, "Separate UMI groups with T bases instead of dashes (e.g. AAAATTAAAATTAAAA)")
	ontPrimersCmd.Flags().BoolVar(&ontTrimFlanking, "trim-flanking", false, "Trim VNP/SSP sequences from passing reads before writing")
	ontPrimersCmd.Flags().BoolVar(&ontSenseCorrect, "sense-correct", false, "Correct the read to be in the sense orientation (SSP+/VNP-)")

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
