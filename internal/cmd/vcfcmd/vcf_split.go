package vcfcmd

import (
	"errors"
	"fmt"
	"io"
	"io/fs"
	"os"
	"strconv"

	"github.com/compgenlab/cghts/vcf"
	"github.com/compgenlab/cgkit/internal/locator"
	"github.com/spf13/cobra"
)

var (
	vcfSplitOut   string
	vcfSplitNum   int
	vcfSplitForce bool
)

// splitChunkPath names the nth chunk of a series. One function so the writer,
// the overwrite guard and the cleanup cannot disagree about the naming, which
// is what vcf-concat --chunks relies on to find them.
func splitChunkPath(base string, n int) string {
	return base + "." + strconv.Itoa(n) + ".vcf.gz"
}

// existingChunks lists the chunks already on disk for a base name, stopping at
// the first gap -- which is exactly how vcf-concat --chunks reads a series, so
// the guard sees the same files the consumer would.
func existingChunks(base string) []string {
	var out []string
	for n := 1; ; n++ {
		p := splitChunkPath(base, n)
		if _, err := os.Stat(p); err != nil {
			return out
		}
		out = append(out, p)
	}
}

var vcfSplitCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-split <input.vcf>",
	Short:       "Split a VCF file into smaller files with N variants each",
	Long: `Split a VCF file into multiple bgzipped files of N variants each. Each output
file gets a fresh copy of the header. Outputs are named BASE.1.vcf.gz,
BASE.2.vcf.gz, and so on. Recombine them with "vcf-concat --chunks BASE.1.vcf.gz".

A series is written or it is not. Re-running over an existing one is refused
without --force, and a run that fails partway removes the chunks it wrote.
Both exist because vcf-concat --chunks stops at the first *missing* file rather
than the first invalid one: a truncated tail chunk, or a leftover chunk from a
longer previous run, is silently reassembled as though it were the whole thing.

  --out BASE    base output name (required)
  --num N       variants per output file (required)
  --force       overwrite an existing series at BASE
  --tbi         also write a tabix index (BASE.N.vcf.gz.tbi) for each chunk`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		if vcfSplitOut == "" {
			return fmt.Errorf("you must specify a base output name with --out")
		}
		if err := locator.CheckLocalOutput("--out", vcfSplitOut); err != nil {
			return err
		}
		if vcfSplitNum <= 0 {
			return fmt.Errorf("--num must be a positive number")
		}

		// Before opening anything: this is the last moment the previous series
		// is still intact.
		if old := existingChunks(vcfSplitOut); len(old) > 0 {
			if !vcfSplitForce {
				return fmt.Errorf("%s already has %d chunk(s), starting with %s; "+
					"pass --force to replace the series", vcfSplitOut, len(old), old[0])
			}
			// Every one of them, not just the ones this run will overwrite. A
			// shorter series would otherwise leave the old tail in place, and
			// vcf-concat would read it as a continuation of the new one.
			for _, p := range old {
				if err := os.Remove(p); err != nil {
					return fmt.Errorf("replacing previous series: %w", err)
				}
				if err := os.Remove(p + ".tbi"); err != nil && !errors.Is(err, fs.ErrNotExist) {
					return fmt.Errorf("replacing previous series: %w", err)
				}
			}
		}

		reader, err := openVcfInput(cmd, args[0])
		if err != nil {
			return err
		}
		defer reader.Close()

		header, err := reader.Header()
		if err != nil {
			return err
		}
		stampVcfProvenance(header, "vcf-split")

		// Each chunk is a separate file, so with --tbi each gets its own index.
		var writer *vcf.VcfWriter
		fileNum, inFile := 0, 0
		chunkPath := ""
		// A refused index does not abort the split. writeVcfTbi refuses to index
		// an unsorted chunk *after* that chunk is complete and valid, and
		// stopping there would leave a partial series -- the very thing
		// vcf-concat --chunks cannot detect. So the series is finished, the
		// indexes are skipped, and the first failure is reported at the end.
		var indexErr error
		closeChunk := func() error {
			if writer == nil {
				return nil
			}
			err := writer.Close()
			writer = nil
			if err != nil || !vcfTbi {
				return err
			}
			if terr := writeVcfTbi(chunkPath); terr != nil && indexErr == nil {
				indexErr = terr
			}
			return nil
		}
		// A failed run leaves nothing. A partial series is the dangerous
		// outcome here, not merely an untidy one -- see the command help.
		done := false
		defer func() {
			if done {
				return
			}
			// nil once a chunk boundary closed it, and VcfWriter.Close does
			// not tolerate a nil receiver.
			if writer != nil {
				writer.Close()
			}
			for n := 1; n <= fileNum; n++ {
				os.Remove(splitChunkPath(vcfSplitOut, n))
				os.Remove(splitChunkPath(vcfSplitOut, n) + ".tbi")
			}
		}()

		for {
			rec, err := reader.NextRecord()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if writer == nil {
				fileNum++
				chunkPath = splitChunkPath(vcfSplitOut, fileNum)
				w, oerr := vcf.OpenVcfWriter(chunkPath)
				if oerr != nil {
					return oerr
				}
				writer = w
				if herr := writer.WriteHeader(header); herr != nil {
					return herr
				}
				inFile = 0
			}
			if err := writer.WriteRecord(rec); err != nil {
				return err
			}
			inFile++
			if inFile >= vcfSplitNum {
				if err := closeChunk(); err != nil {
					return err
				}
			}
		}
		if err := closeChunk(); err != nil {
			return err
		}
		// The series is whole, so it is kept even when an index was refused.
		done = true
		return indexErr
	},
}

func init() {
	f := vcfSplitCmd.Flags()
	f.StringVar(&vcfSplitOut, "out", "", "Base output name (outputs are BASE.N.vcf.gz)")
	f.IntVar(&vcfSplitNum, "num", 0, "Number of variants per output file")
	f.BoolVar(&vcfSplitForce, "force", false, "Overwrite an existing series at --out")
	f.BoolVar(&vcfTbi, "tbi", false, "Also write a tabix (.tbi) index for each chunk")
}
