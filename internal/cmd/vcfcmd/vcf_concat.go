package vcfcmd

import (
	"fmt"
	"io"

	"github.com/compgenlab/hts/vcf"
	"github.com/spf13/cobra"
)

var (
	vcfConcatOutput string
	vcfConcatChunks bool
)

var vcfConcatCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": sinceVersion},
	Use:         "vcf-concat <input1.vcf> [input2.vcf ...]",
	Short:       "Concatenate VCF files with the same samples but different variants",
	Long: `Concatenate coordinate-sorted VCF files that share the same samples into one
sorted stream. Records are merged by contig order (from the ##contig header lines)
and position; overlapping positions are an error. Header INFO/FORMAT/FILTER/ALT
definitions from all inputs are combined, and every input's chromosomes must be
declared as ##contig lines in a consistent order.

  --chunks    treat each argument as the first file of a numbered chunk sequence
              (base.1.vcf.gz, base.2.vcf.gz, ...) and read it one file at a time,
              so recombining thousands of vcf-split chunks stays within the
              file-descriptor limit`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}

		var streams []*recordSource
		closeAll := func() {
			for _, s := range streams {
				if s != nil {
					s.close()
				}
			}
		}
		for _, fn := range args {
			s, err := openConcatStream(cmd, fn, vcfConcatChunks)
			if err != nil {
				closeAll()
				return err
			}
			streams = append(streams, s)
		}
		defer closeAll()

		header := streams[0].header
		for _, s := range streams[1:] {
			if err := unionHeaderInto(header, s.header, false); err != nil {
				return err
			}
		}
		stampVcfProvenance(header, "vcf-concat")

		order := map[string]int{}
		for i, c := range header.ContigNames() {
			order[c] = i
		}

		writer, closeFn, err := openVcfWriter(cmd, vcfConcatOutput)
		if err != nil {
			return err
		}
		if err := writer.WriteHeader(header); err != nil {
			return err
		}

		cur := make([]*vcf.VcfRecord, len(streams))
		exhausted := make([]bool, len(streams))
		load := func(i int) error {
			if cur[i] != nil || exhausted[i] {
				return nil
			}
			rec, err := streams[i].next()
			if err == io.EOF {
				exhausted[i] = true
				return nil
			}
			if err != nil {
				return err
			}
			if _, ok := order[rec.Chrom]; !ok {
				return fmt.Errorf("unknown chromosome: %s (all chromosomes must be listed as a ##contig and in the same order across inputs)", rec.Chrom)
			}
			cur[i] = rec
			return nil
		}

		for {
			for i := range streams {
				if err := load(i); err != nil {
					return err
				}
			}
			low := -1
			for i := range streams {
				if cur[i] == nil {
					continue
				}
				if low == -1 {
					low = i
					continue
				}
				ci, cl := order[cur[i].Chrom], order[cur[low].Chrom]
				switch {
				case ci < cl, ci == cl && cur[i].Pos < cur[low].Pos:
					low = i
				case ci == cl && cur[i].Pos == cur[low].Pos:
					return fmt.Errorf("overlapping variant positions found: %s:%d", cur[i].Chrom, cur[i].Pos)
				}
			}
			if low == -1 {
				break
			}
			if err := writer.WriteRecord(cur[low]); err != nil {
				return err
			}
			cur[low] = nil
		}

		if closeFn != nil {
			return closeFn()
		}
		return writer.Close()
	},
}

// openConcatStream opens filename as a record source. With chunks set, it reads
// the numbered chunk sequence beginning at filename one file at a time.
func openConcatStream(cmd *cobra.Command, filename string, chunks bool) (*recordSource, error) {
	if chunks {
		c, err := vcf.NewChunkedVcfReader(filename)
		if err != nil {
			return nil, err
		}
		header, err := c.Header()
		if err != nil {
			c.Close()
			return nil, err
		}
		return &recordSource{header: header, next: c.NextRecord, close: c.Close}, nil
	}
	reader, err := openVcfInput(cmd, filename)
	if err != nil {
		return nil, err
	}
	header, err := reader.Header()
	if err != nil {
		reader.Close()
		return nil, err
	}
	return &recordSource{
		header: header,
		next:   reader.NextRecord,
		close:  func() error { reader.Close(); return nil },
	}, nil
}

// unionHeaderInto folds src's definitions into dst (first file wins). Samples in
// src must already exist in dst. When addContigs is false (vcf-concat) a contig
// in src that is absent from dst is an error; when true (vcf-merge) it is added.
func unionHeaderInto(dst, src *vcf.VcfHeader, addContigs bool) error {
	dstSamples := map[string]bool{}
	for _, s := range dst.Samples() {
		dstSamples[s] = true
	}
	for _, s := range src.Samples() {
		if !dstSamples[s] {
			return fmt.Errorf("file contains an extra sample: %s", s)
		}
	}
	for _, id := range src.FilterIDs() {
		if _, ok := dst.FilterDef(id); !ok {
			d, _ := src.FilterDef(id)
			dst.AddFilter(d)
		}
	}
	for _, id := range src.InfoIDs() {
		if _, ok := dst.InfoDef(id); !ok {
			d, _ := src.InfoDef(id)
			dst.AddInfo(d)
		}
	}
	for _, id := range src.FormatIDs() {
		if _, ok := dst.FormatDef(id); !ok {
			d, _ := src.FormatDef(id)
			dst.AddFormat(d)
		}
	}
	for _, id := range src.ContigNames() {
		if _, ok := dst.ContigDef(id); !ok {
			if !addContigs {
				return fmt.Errorf("unknown contig: %s (all chromosomes must be listed as a ##contig and in the same order across inputs)", id)
			}
			d, _ := src.ContigDef(id)
			dst.AddContig(d)
		}
	}
	for _, id := range src.AltIDs() {
		if _, ok := dst.AltDef(id); !ok {
			d, _ := src.AltDef(id)
			dst.AddAlt(d)
		}
	}
	existing := map[string]bool{}
	for _, l := range dst.OtherLines() {
		existing[l] = true
	}
	for _, l := range src.OtherLines() {
		if !existing[l] {
			dst.AddLine(l)
			existing[l] = true
		}
	}
	return nil
}

func init() {
	f := vcfConcatCmd.Flags()
	f.StringVarP(&vcfConcatOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
	f.BoolVar(&vcfConcatChunks, "chunks", false, "Treat each argument as the first file of a numbered chunk sequence (base.1.vcf.gz, ...)")
}
