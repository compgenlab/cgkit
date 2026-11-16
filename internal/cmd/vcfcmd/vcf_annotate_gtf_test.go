package vcfcmd

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio/tabix"
)

// disagreeingGTF is one + strand gene G whose transcripts disagree about
// position 151: it is in the non-coding transcript's exon, and upstream of the
// coding transcript's CDS without being in one of its exons. nc and cd name
// the two transcripts. It has gene and transcript rows, as GENCODE and Ensembl
// do: the indexed reader takes a gene's extent from them, and without them can
// see only the rows overlapping the position.
func disagreeingGTF(nc, cd string) []string {
	g := `gene_id "G"; gene_name "GENE"; transcript_id "`
	return []string{
		"chr1\tt\tgene\t101\t600\t.\t+\t.\tgene_id \"G\"; gene_name \"GENE\";",
		"chr1\tt\ttranscript\t101\t200\t.\t+\t.\t" + g + nc + `";`,
		"chr1\tt\ttranscript\t301\t600\t.\t+\t.\t" + g + cd + `";`,
		"chr1\tt\texon\t101\t200\t.\t+\t.\t" + g + nc + `";`,
		"chr1\tt\texon\t301\t400\t.\t+\t.\t" + g + cd + `";`,
		"chr1\tt\tCDS\t321\t400\t.\t+\t0\t" + g + cd + `";`,
		"chr1\tt\texon\t501\t600\t.\t+\t.\t" + g + cd + `";`,
		"chr1\tt\tCDS\t501\t550\t.\t+\t0\t" + g + cd + `";`,
	}
}

// GTF_REGION used to depend on the transcripts' names: cghts visited them in
// sorted-ID order and carried flags from one to the next, so named coding-first
// the position came out a non-coding exon and named the other way a 5'UTR
// (compgenlab/cghts#61). It is exonic in one transcript and before the other's
// CDS, so it is a 5'UTR, whatever the names and whether the GTF is read whole
// or through its index.
func TestVcfAnnotateGtfRegionIgnoresTranscriptNames(t *testing.T) {
	vcf := filepath.Join(t.TempDir(), "v.vcf")
	if err := os.WriteFile(vcf, []byte("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t151\t.\tA\tG\t.\tPASS\t.\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	for _, names := range [][2]string{{"Z_nc", "A_cd"}, {"A_nc", "Z_cd"}} {
		lines := disagreeingGTF(names[0], names[1])
		dir := t.TempDir()

		plain := filepath.Join(dir, "plain.gtf")
		if err := os.WriteFile(plain, []byte(strings.Join(lines, "\n")+"\n"), 0o644); err != nil {
			t.Fatal(err)
		}
		indexed := filepath.Join(dir, "indexed.gtf.gz")
		w := tabix.NewWriter(indexed, tabix.NewWriterOpts().GFF().AutoIndex())
		for _, l := range lines {
			if err := w.Write(l); err != nil {
				t.Fatal(err)
			}
		}
		if err := w.Close(); err != nil {
			t.Fatal(err)
		}

		for _, gtf := range []string{plain, indexed} {
			out := runVcf(t, "vcf-annotate", "--gtf", gtf, vcf)
			var row string
			for _, l := range strings.Split(out, "\n") {
				if strings.HasPrefix(l, "chr1\t151\t") {
					row = l
				}
			}
			if !strings.Contains(row, "GTF_REGION=5_utr") || !strings.Contains(row, "GTF_GENE=GENE") {
				t.Errorf("transcripts %v, %s: want GTF_GENE=GENE and GTF_REGION=5_utr, got\n%s", names, filepath.Base(gtf), row)
			}
		}
	}
}
