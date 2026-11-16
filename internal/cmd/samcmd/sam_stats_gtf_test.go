package samcmd

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio"
	"github.com/compgenlab/cghts/htsio/tabix"
)

// gtfFixture is one + strand coding gene, 0-based:
//
//	exon1 [1000,1200)  5'UTR [1000,1100) + CDS [1100,1200)
//	intron [1200,1500)  coding intron (between the CDS bounds)
//	exon2 [1500,2000)  CDS [1500,1800) + 3'UTR [1800,2000)
//
// It has one transcript, so nothing in it depends on how transcripts that
// disagree are reconciled -- that is cghts's to test.
const gtfFixture = "" +
	"chr1\tt\tgene\t1001\t2000\t.\t+\t.\tgene_id \"G\"; gene_name \"GENE\"; gene_type \"protein_coding\";\n" +
	"chr1\tt\ttranscript\t1001\t2000\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";\n" +
	"chr1\tt\texon\t1001\t1200\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";\n" +
	"chr1\tt\tCDS\t1101\t1200\t.\t+\t0\tgene_id \"G\"; transcript_id \"T\";\n" +
	"chr1\tt\texon\t1501\t2000\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";\n" +
	"chr1\tt\tCDS\t1501\t1800\t.\t+\t0\tgene_id \"G\"; transcript_id \"T\";\n"

// gtfReads lands one read in each region worth telling apart. A read is
// classified at its 5' end: a forward read's first base, a reverse read's
// alignment end (ngsutilsj's position, one past its last base).
var gtfReads = []string{
	"coding\t0\tchr1\t1110\t60\t50M",          // CDS, sense
	"utr5\t0\tchr1\t1010\t60\t50M",            // before the CDS
	"cintron\t0\tchr1\t1300\t60\t50M",         // between the CDS bounds
	"utr3\t0\tchr1\t1900\t60\t50M",            // after the CDS
	"intergenic\t0\tchr1\t50000\t60\t50M",     // no gene
	"anti\t16\tchr1\t1121\t60\t50M",           // reverse: 5' end in the CDS, wrong strand
	"junction\t0\tchr1\t1150\t60\t25M300N25M", // spliced, starting in the CDS
	"mito\t0\tchrM\t100\t60\t50M",             // chrM is its own category
}

func writeGtfFixture(t *testing.T) (gtfPath, samPath string) {
	t.Helper()
	dir := t.TempDir()
	gtfPath = filepath.Join(dir, "genes.gtf.gz")
	w := tabix.NewWriter(gtfPath, tabix.NewWriterOpts().GFF().AutoIndex())
	for _, line := range strings.Split(strings.TrimSuffix(gtfFixture, "\n"), "\n") {
		if err := w.Write(line); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	var b strings.Builder
	b.WriteString("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chrM\tLN:16569\n")
	for _, r := range gtfReads {
		b.WriteString(r + "\t*\t0\t0\t" + strings.Repeat("A", 50) + "\t" + strings.Repeat("I", 50) + "\tNH:i:1\n")
	}
	samPath = filepath.Join(dir, "reads.sam")
	if err := os.WriteFile(samPath, []byte(b.String()), 0o644); err != nil {
		t.Fatal(err)
	}
	return gtfPath, samPath
}

func TestSamStatsGtf(t *testing.T) {
	gtfPath, samPath := writeGtfFixture(t)
	genes, closeGenes, err := openGeneModel(gtfPath)
	if err != nil {
		t.Fatal(err)
	}
	defer closeGenes()
	reader, err := htsio.NewSamReader(samPath)
	if err != nil {
		t.Fatal(err)
	}
	defer reader.Close()
	header, _ := reader.Header()
	opt := statsOpts{genes: genes}
	res, err := computeStats(reader, header, opt, io.Discard)
	if err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	writeStats(&out, res, opt)
	got := out.String()
	got = got[strings.Index(got, "[Gene regions]"):]
	if got != wantGeneSection {
		t.Errorf("gene section:\n%s\nwant:\n%s", got, wantGeneSection)
	}
}

// wantGeneSection is ngsutilsj bam-stats --gtf's own output for this fixture,
// byte for byte.
const wantGeneSection = "[Gene regions]\n" +
	"Junction                        \t1\n" +
	"Coding                          \t1\n" +
	"5'UTR                           \t1\n" +
	"3'UTR                           \t1\n" +
	"Non-coding junction             \t0\n" +
	"Non-coding exon                 \t0\n" +
	"Coding intron                   \t1\n" +
	"5'UTR intron                    \t0\n" +
	"3'UTR intron                    \t0\n" +
	"Non-coding intron               \t0\n" +
	"Intergenic                      \t1\n" +
	"Mitochondrial                   \t1\n" +
	"Junction (anti-sense)           \t0\n" +
	"Coding (anti-sense)             \t1\n" +
	"5'UTR (anti-sense)              \t0\n" +
	"3'UTR (anti-sense)              \t0\n" +
	"Non-coding junction (anti-sense)\t0\n" +
	"Non-coding exon (anti-sense)    \t0\n" +
	"Coding intron (anti-sense)      \t0\n" +
	"5'UTR intron (anti-sense)       \t0\n" +
	"3'UTR intron (anti-sense)       \t0\n" +
	"Non-coding intron (anti-sense)  \t0\n" +
	"\n" +
	"Orientation\tunstranded\n" +
	"Sense      \t5\n" +
	"Anti-sense \t1\n" +
	"\n" +
	"Exonic     \t5\n" +
	"Intronic   \t1\n" +
	"\n" +
	"Junction   \t1\n" +
	"\n" +
	"Sense/anti-sense ratio     \t5.00\n" +
	"Exonic/intronic ratio      \t5.00\n" +
	"Exonic/genomic ratio       \t5.00\n" +
	"Non-junction/junction ratio\t4.00\n"

// A GTF with no index is refused rather than loaded whole.
func TestSamStatsGtfNeedsAnIndex(t *testing.T) {
	path := filepath.Join(t.TempDir(), "genes.gtf")
	if err := os.WriteFile(path, []byte(gtfFixture), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, _, err := openGeneModel(path); err == nil || !strings.Contains(err.Error(), "no tabix index") {
		t.Errorf("err = %v, want a missing-index refusal", err)
	}
}
