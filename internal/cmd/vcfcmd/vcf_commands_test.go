package vcfcmd

import (
	"bytes"
	"strings"
	"testing"
	"time"

	"github.com/compgenlab/cgio/internal/buildinfo"
	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// runVcf executes a vcf subcommand against a fresh root and returns its stdout.
// Command flags are bound to package globals, so they are reset to their
// defaults before each run to keep tests independent of ordering.
func runVcf(t *testing.T, args ...string) string {
	t.Helper()
	root := &cobra.Command{Use: "cgio"}
	InitCmd(root)
	for _, c := range root.Commands() {
		c.Flags().VisitAll(func(f *pflag.Flag) {
			_ = f.Value.Set(f.DefValue)
			f.Changed = false
		})
	}
	// Array-backed flags append on Set, so the default reset above does not
	// clear them; reset their globals directly.
	vcfExportInfo = nil
	vcfExportFormat = nil
	vcfReorderSamples = nil
	vcfReorderSamplesFile = ""
	vcfStatsInfoTally = nil
	vcfStatsInfoPresent = nil
	vcfAnnotateTags = nil
	vcfAnnotateBed = nil
	vcfAnnotateBedFlag = nil
	vcfAnnotateFormatBed = nil
	vcfAnnotateTab = nil
	vcfAnnotateFormatTab = nil
	vcfAnnotateVcf = nil
	vcfAnnotateVcfFlag = nil
	vcfAnnotateInFile = nil
	var buf bytes.Buffer
	root.SetOut(&buf)
	root.SetErr(&buf)
	root.SetArgs(args)
	if err := root.Execute(); err != nil {
		t.Fatalf("Execute(%v): %v", args, err)
	}
	return buf.String()
}

func TestVcfSamples(t *testing.T) {
	want := "NORMAL\nTUMOR\n"
	if got := runVcf(t, "vcf-samples", "testdata/sample.vcf"); got != want {
		t.Errorf("vcf-samples mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestVcfToBed(t *testing.T) {
	want := "chr1\t99\t100\tSNV\n" +
		"chr1\t199\t200\tSNV\n" +
		"chr1\t299\t300\tINS\n" +
		"chr2\t499\t900\tDEL\n"
	if got := runVcf(t, "vcf-tobed", "testdata/sample.vcf"); got != want {
		t.Errorf("vcf-tobed mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestVcfToBedPassing(t *testing.T) {
	want := "chr1\t99\t100\tSNV\n" +
		"chr1\t299\t300\tINS\n" +
		"chr2\t499\t900\tDEL\n"
	if got := runVcf(t, "vcf-tobed", "--passing", "testdata/sample.vcf"); got != want {
		t.Errorf("vcf-tobed --passing mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestVcfExport(t *testing.T) {
	want := "chrom\tpos\tref\talt\tID\tQUAL\tFILTER\tDP\tAF\tNORMAL:AD\tTUMOR:AD\n" +
		"chr1\t100\tA\tG\trs1\t50.0\tPASS\t30\t0.5\t28,2\t15,15\n" +
		"chr1\t200\tC\tT\t.\t20.0\tlowqual\t10\t0.1\t9,1\t5,5\n" +
		"chr1\t300\tG\tGA\t.\t40.0\tPASS\t25\t\t25,0\t12,13\n" +
		"chr2\t500\tA\t<DEL>\t.\t99.0\tPASS\t\t\t\t\n" +
		"chr2\t1000\tT\tT[chr5:2000[\tbnd1\t60.0\tPASS\t\t\t\t\n"
	got := runVcf(t, "vcf-export", "--no-vcf-header",
		"--id", "--qual", "--filter", "--info", "DP", "--info", "AF", "--format", "AD",
		"testdata/sample.vcf")
	if got != want {
		t.Errorf("vcf-export mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestVcfExportAlleleSelectors(t *testing.T) {
	want := "chrom\tpos\tref\talt\tTUMOR:AD\tTUMOR:AD\tTUMOR:AD\n" +
		"chr1\t100\tA\tG\t30\t15\t15\n" +
		"chr1\t200\tC\tT\t10\t5\t5\n" +
		"chr1\t300\tG\tGA\t25\t12\t13\n" +
		"chr2\t500\tA\t<DEL>\t\t\t\n" +
		"chr2\t1000\tT\tT[chr5:2000[\t\t\t\n"
	got := runVcf(t, "vcf-export", "--no-vcf-header",
		"--format", "AD:TUMOR:sum", "--format", "AD:TUMOR:ref", "--format", "AD:TUMOR:alt1",
		"testdata/sample.vcf")
	if got != want {
		t.Errorf("vcf-export allele selectors mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestVcfExportOnlySnvsMissingBlank(t *testing.T) {
	want := "chrom\tpos\tref\talt\tAF\n" +
		"chr1\t100\tA\tG\t0.5\n" +
		"chr1\t200\tC\tT\t0.1\n"
	got := runVcf(t, "vcf-export", "--no-vcf-header", "--only-snvs", "--missing-blank",
		"--info", "AF", "testdata/sample.vcf")
	if got != want {
		t.Errorf("vcf-export --only-snvs --missing-blank mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestVcfToBedRegion(t *testing.T) {
	// 1-based inclusive region over the tabix-indexed file.
	want := "chr1\t199\t200\tSNV\n" +
		"chr1\t299\t300\tINS\n"
	if got := runVcf(t, "vcf-tobed", "--region", "chr1:200-1000", "testdata/sample.vcf.gz"); got != want {
		t.Errorf("vcf-tobed --region mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestVcfToBedRegionWholeChrom(t *testing.T) {
	want := "chr2\t499\t900\tDEL\n"
	if got := runVcf(t, "vcf-tobed", "--region", "chr2", "testdata/sample.vcf.gz"); got != want {
		t.Errorf("vcf-tobed --region chr2 mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestVcfExportRegion(t *testing.T) {
	want := "chrom\tpos\tref\talt\tDP\n" +
		"chr1\t300\tG\tGA\t25\n"
	got := runVcf(t, "vcf-export", "--no-vcf-header", "--region", "chr1:300-300",
		"--info", "DP", "testdata/sample.vcf.gz")
	if got != want {
		t.Errorf("vcf-export --region mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestVcfRegionRequiresFile(t *testing.T) {
	root := &cobra.Command{Use: "cgio"}
	InitCmd(root)
	var buf bytes.Buffer
	root.SetOut(&buf)
	root.SetErr(&buf)
	root.SetArgs([]string{"vcf-tobed", "--region", "chr1", "-"})
	if err := root.Execute(); err == nil {
		t.Fatal("expected an error when --region is used with stdin")
	}
}

func TestVcfReorder(t *testing.T) {
	// Samples swapped in the #CHROM line and in each record; non-sample columns
	// (including QUAL "50.0") are preserved verbatim.
	got := runVcf(t, "vcf-reorder", "-s", "TUMOR", "-s", "NORMAL", "testdata/sample.vcf")
	wantChrom := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n"
	if !strings.Contains(got, wantChrom) {
		t.Errorf("vcf-reorder #CHROM line missing.\n got: %q", got)
	}
	wantRow := "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30;AF=0.5;DB\tGT:AD\t0/1:15,15\t0/0:28,2\n"
	if !strings.Contains(got, wantRow) {
		t.Errorf("vcf-reorder first data row mismatch.\n got: %q", got)
	}
}

func TestVcfReorderProvenance(t *testing.T) {
	buildinfo.Now = func() time.Time { return time.Date(2020, 1, 2, 0, 0, 0, 0, time.UTC) }
	defer func() { buildinfo.Now = time.Now }()

	got := runVcf(t, "vcf-reorder", "-s", "NORMAL", "-s", "TUMOR", "testdata/sample.vcf")
	if !strings.Contains(got, "##fileDate=20200102\n") {
		t.Errorf("vcf-reorder missing/incorrect ##fileDate.\n got: %q", got)
	}
	if !strings.Contains(got, "##cgio_vcf-reorderCommand=") {
		t.Errorf("vcf-reorder missing provenance command line.\n got: %q", got)
	}
	if !strings.Contains(got, "##cgio_vcf-reorderVersion=dev\n") {
		t.Errorf("vcf-reorder missing provenance version.\n got: %q", got)
	}
}

func TestVcfReorderSubsetByNumber(t *testing.T) {
	// Select only sample 2 (TUMOR) by 1-based number.
	got := runVcf(t, "vcf-reorder", "-s", "2", "testdata/sample.vcf")
	wantChrom := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t2\n"
	if !strings.Contains(got, wantChrom) {
		t.Errorf("vcf-reorder subset #CHROM line mismatch.\n got: %q", got)
	}
	wantRow := "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30;AF=0.5;DB\tGT:AD\t0/1:15,15\n"
	if !strings.Contains(got, wantRow) {
		t.Errorf("vcf-reorder subset data row mismatch.\n got: %q", got)
	}
}

func TestVcfAnnotateGroupA(t *testing.T) {
	out := runVcf(t, "vcf-annotate", "--indel", "--tstv", "--dosage", "--auto-id",
		"testdata/annotate.vcf")
	for _, want := range []string{
		"chr1\t100\tchr1_100_A_G\tA\tG\t50\tPASS\tDP=30;CG_TSTV=TS\tGT:AD:SAC:CG_DS\t0/0:14,1:15,13,1,1:0\t0/1:15,15:8,7,8,7:1\n",
		"chr1\t150\tchr1_150_A_C\tA\tC\t40\tPASS\tDP=25;CG_TSTV=TV\tGT:AD:SAC:CG_DS\t0/0:25,0:13,12,0,0:0\t1/1:0,30:0,0,15,15:2\n",
		"chr1\t300\tchr1_300_G_GA\tG\tGA\t40\tPASS\tDP=25;CG_INSERT;CG_INSLEN=1;CG_INDELLEN=1\tGT:AD:SAC:CG_DS\t0/0:25,0:13,12,0,0:0\t0/1:12,13:6,6,7,6:1\n",
		"chr2\t500\tchr2_500_CAT_C\tCAT\tC\t99\tPASS\tDP=40;CG_DELETE;CG_DELLEN=2;CG_INDELLEN=-2\tGT:AD:SAC:CG_DS\t0/0:40,0:20,20,0,0:0\t0/1:20,20:10,10,11,9:1\n",
	} {
		if !strings.Contains(out, want) {
			t.Errorf("vcf-annotate group A missing line:\n%q\nfull output:\n%s", want, out)
		}
	}
	// New INFO/FORMAT defs declared in the header.
	for _, def := range []string{
		"##INFO=<ID=CG_TSTV,",
		"##INFO=<ID=CG_INSERT,",
		"##FORMAT=<ID=CG_DS,",
	} {
		if !strings.Contains(out, def) {
			t.Errorf("vcf-annotate header missing def %q", def)
		}
	}
}

func TestVcfAnnotateBedTab(t *testing.T) {
	out := runVcf(t, "vcf-annotate",
		"--bed", "REGION:testdata/regions.bed.gz",
		"--bed-flag", "INREG:testdata/regions.bed.gz",
		"--tab", "SCORE:testdata/scores.tab.gz,5,n",
		"--tab", "LBL:testdata/scores.tab.gz,6,alt=4",
		"testdata/annotate.vcf")
	// chr1:100 A>G overlaps promoterX, score row A>G (0.91, hot).
	if !strings.Contains(out, "DP=30;REGION=promoterX;INREG;SCORE=0.91;LBL=hot\t") {
		t.Errorf("chr1:100 annotations wrong:\n%s", out)
	}
	// chr1:150 A>C: region promoterX, score row A>C (0.40, cold).
	if !strings.Contains(out, "DP=25;REGION=promoterX;INREG;SCORE=0.40;LBL=cold\t") {
		t.Errorf("chr1:150 annotations wrong:\n%s", out)
	}
	// chr1:300 indel in exonY, no score row there.
	if !strings.Contains(out, "DP=25;REGION=exonY;INREG\t") {
		t.Errorf("chr1:300 annotations wrong:\n%s", out)
	}
	// Header defs added (tabix-style).
	for _, def := range []string{"##INFO=<ID=REGION,", "##INFO=<ID=INREG,", "##INFO=<ID=SCORE,", "##INFO=<ID=LBL,"} {
		if !strings.Contains(out, def) {
			t.Errorf("missing header def %q", def)
		}
	}
}

func TestVcfAnnotateVcfSource(t *testing.T) {
	out := runVcf(t, "vcf-annotate",
		"--vcf", "KAF:AF:testdata/source.vcf.gz",
		"--vcf-id", "testdata/source.vcf.gz",
		"--vcf-flag", "KNOWN:testdata/source.vcf.gz",
		"testdata/annotate.vcf")
	// chr1:100 A>G matches rsA: ID copied, KAF pulled, KNOWN flag.
	if !strings.Contains(out, "chr1\t100\trsA\tA\tG\t50\tPASS\tDP=30;KAF=0.20;KNOWN\t") {
		t.Errorf("chr1:100 vcf annotations wrong:\n%s", out)
	}
	// chr2:500 deletion matches rsC (exact ref/alt): ID + flag, no AF in source.
	if !strings.Contains(out, "chr2\t500\trsC\tCAT\tC\t99\tPASS\tDP=40;KNOWN\t") {
		t.Errorf("chr2:500 vcf annotations wrong:\n%s", out)
	}
	for _, def := range []string{"##INFO=<ID=KAF,", "##INFO=<ID=KNOWN,"} {
		if !strings.Contains(out, def) {
			t.Errorf("missing header def %q", def)
		}
	}
}

func TestVcfAnnotateVcfPassing(t *testing.T) {
	// rsB (chr1:150) is "rej"; with @ it is skipped.
	out := runVcf(t, "vcf-annotate", "--vcf", "KAF:AF:testdata/source.vcf.gz:@", "testdata/annotate.vcf")
	if strings.Contains(out, "KAF=0.10") {
		t.Errorf("passing-only should skip the rejected rsB record:\n%s", out)
	}
	if !strings.Contains(out, "KAF=0.20") {
		t.Errorf("passing-only dropped a passing record:\n%s", out)
	}
}

func TestVcfAnnotateInFile(t *testing.T) {
	// Flag mode: DP values 30 and 25 are in dp_set.txt.
	out := runVcf(t, "vcf-annotate", "--in-file", "HIT:DP:testdata/dp_set.txt", "testdata/annotate.vcf")
	if !strings.Contains(out, "DP=30;HIT\t") || !strings.Contains(out, "DP=25;HIT\t") {
		t.Errorf("in-file flag mismatch:\n%s", out)
	}
	if strings.Contains(out, "DP=10;HIT") || strings.Contains(out, "DP=40;HIT") {
		t.Errorf("in-file flagged a non-listed value:\n%s", out)
	}
	// tabcol mode: value from column 2.
	out = runVcf(t, "vcf-annotate", "--in-file", "LABEL:DP:testdata/dp_val.txt:tabcol=2", "testdata/annotate.vcf")
	if !strings.Contains(out, "DP=30;LABEL=common\t") || !strings.Contains(out, "DP=25;LABEL=rare\t") {
		t.Errorf("in-file tabcol mismatch:\n%s", out)
	}
}

func TestVcfAnnotateColumnByName(t *testing.T) {
	// scores_hdr.tab.gz has a header; columns referenced by name. Exact ref+alt
	// matching also reaches the chr2:500 deletion.
	out := runVcf(t, "vcf-annotate",
		"--tab", "LBL:testdata/scores_hdr.tab.gz,label,alt=alt,ref=ref",
		"testdata/annotate.vcf")
	if !strings.Contains(out, "DP=30;LBL=hot\t") {
		t.Errorf("chr1:100 LBL missing:\n%s", out)
	}
	if !strings.Contains(out, "DP=40;LBL=sv\t") {
		t.Errorf("chr2:500 (deletion, ref+alt match) LBL missing:\n%s", out)
	}
	if !strings.Contains(out, "##INFO=<ID=LBL,") {
		t.Errorf("LBL header def missing")
	}
}

func TestVcfAnnotateFormatBed(t *testing.T) {
	// Annotate sample 0 (NORMAL): the FORMAT keys are derived from sample 0, so
	// the new field surfaces; TUMOR lacks it and its trailing value is trimmed.
	out := runVcf(t, "vcf-annotate",
		"--format-bed", "REGION:NORMAL:testdata/regions.bed.gz",
		"testdata/annotate.vcf")
	if !strings.Contains(out, "GT:AD:SAC:REGION\t0/0:14,1:15,13,1,1:promoterX\t0/1:15,15:8,7,8,7\n") {
		t.Errorf("format-bed output wrong:\n%s", out)
	}
}

func TestVcfAnnotateProvenance(t *testing.T) {
	buildinfo.Now = func() time.Time { return time.Date(2021, 3, 4, 0, 0, 0, 0, time.UTC) }
	defer func() { buildinfo.Now = time.Now }()
	out := runVcf(t, "vcf-annotate", "--tstv", "testdata/annotate.vcf")
	if !strings.Contains(out, "##fileDate=20210304\n") {
		t.Errorf("vcf-annotate missing ##fileDate")
	}
	if !strings.Contains(out, "##cgio_vcf-annotateCommand=") {
		t.Errorf("vcf-annotate missing provenance command")
	}
}

func TestVcfAnnotatePassing(t *testing.T) {
	// --passing drops the lowqual variant at chr1:140.
	out := runVcf(t, "vcf-annotate", "--tstv", "--passing", "testdata/annotate.vcf")
	if strings.Contains(out, "\t140\t") {
		t.Errorf("vcf-annotate --passing should drop chr1:140\n%s", out)
	}
	if !strings.Contains(out, "\t100\t") {
		t.Errorf("vcf-annotate --passing dropped a passing variant\n%s", out)
	}
}

func TestVcfStats(t *testing.T) {
	want := "Total variants:\t5\n" +
		"Filtered variants:\t1\n" +
		"Passing variants:\t4\n" +
		"\n" +
		"SNV:\t2\n" +
		"Indels:\t3\n" +
		"Reference-only:\t0\n" +
		"\n" +
		"Transitions:\t2\n" +
		"Transversions:\t0\n" +
		"Ts/Tv ratio:\t\n" +
		"\n" +
		"[Filters]\n" +
		"lowqual: 1\n"
	if got := runVcf(t, "vcf-stats", "testdata/sample.vcf"); got != want {
		t.Errorf("vcf-stats mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestVcfStatsInfoTally(t *testing.T) {
	want := "[SVTYPE]\n" +
		"BND\t1\n" +
		"DEL\t1\n" +
		"*missing*\t3\n" +
		"\n" +
		"[DB]\n" +
		"Present\t1\n" +
		"Absent\t4\n"
	got := runVcf(t, "vcf-stats", "--info-tally", "SVTYPE", "--info-present", "DB", "testdata/sample.vcf")
	if !strings.HasSuffix(got, want) {
		t.Errorf("vcf-stats info-tally mismatch.\n got: %q\nwant suffix: %q", got, want)
	}
}

func TestVcfTsTv(t *testing.T) {
	want := "Transitions (Ts)\t2\n" +
		"Transversions (Tv)\t0\n" +
		"Ts/Tv ratio\tInfinity\n"
	if got := runVcf(t, "vcf-tstv", "testdata/sample.vcf"); got != want {
		t.Errorf("vcf-tstv mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestVcfToBedIncludePosPadding(t *testing.T) {
	want := "chr1\t94\t105\tchr1_100\n" +
		"chr1\t194\t205\tchr1_200\n" +
		"chr1\t294\t305\tchr1_300\n" +
		"chr2\t494\t905\tchr2_500\n"
	if got := runVcf(t, "vcf-tobed", "--include-pos", "--padding", "5", "testdata/sample.vcf"); got != want {
		t.Errorf("vcf-tobed --include-pos --padding mismatch.\n got: %q\nwant: %q", got, want)
	}
}
