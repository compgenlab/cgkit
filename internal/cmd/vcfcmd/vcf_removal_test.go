package vcfcmd

import (
	"strings"
	"testing"
)

// vcf-chrfix and vcf-strip both drop a *set* of header entries by ranging over
// an accessor and removing as they go. That used to be unsafe: cghts compacted
// the order slice in place, and since range fixes the length up front while the
// elements slide down past the cursor, every second entry was skipped. Both
// commands worked around it by copying the accessor's result first.
//
// cghts v0.10.3 makes the removers allocate, so the copies are gone. These tests
// are what makes that safe rather than merely plausible: each removes a *run* of
// consecutive entries, which is the shape that exposes the skip. Neither command
// had a test over this path before.
//
// testdata/manycontigs.vcf exists for exactly this -- five consecutive
// non-primary contigs, and four consecutive INFO/FORMAT/FILTER definitions, so a
// stride-2 bug leaves survivors rather than happening to come out right.

func TestChrFixDropsEveryNonPrimaryContig(t *testing.T) {
	out := runVcf(t, "vcf-chrfix", "--primary-human", "testdata/manycontigs.vcf")

	for _, gone := range []string{"chrUn_a", "chrUn_b", "chrUn_c", "chrUn_d", "chrUn_e"} {
		if strings.Contains(out, "ID="+gone) {
			t.Errorf("contig %s survived --primary-human:\n%s", gone, headerOf(out))
		}
		// The records for those contigs go too, not just the declarations.
		for _, line := range strings.Split(out, "\n") {
			if strings.HasPrefix(line, gone+"\t") {
				t.Errorf("a record on dropped contig %s survived: %q", gone, line)
			}
		}
	}
	// The primary contigs stay, or the test proves only that everything went.
	for _, kept := range []string{"chr1", "chr2"} {
		if !strings.Contains(out, "##contig=<ID="+kept+",") {
			t.Errorf("primary contig %s was dropped:\n%s", kept, headerOf(out))
		}
	}
}

// The rename path removes and re-adds inside the same loop, so it both shrinks
// and grows the order slice while iterating it -- a harder case than plain
// removal. Every one of the seven contigs is renamed here, so a skipped entry
// leaves a UCSC spelling behind in an otherwise-Ensembl header.
func TestChrFixRenamesEveryContig(t *testing.T) {
	out := runVcf(t, "vcf-chrfix", "--ensembl", "testdata/manycontigs.vcf")
	header := headerOf(out)

	for _, tc := range []struct{ old, want string }{
		{"chr1", "1"},
		{"chr2", "2"},
		// --ensembl strips the prefix from every contig, not only the ones with
		// a known Ensembl equivalent.
		{"chrUn_a", "Un_a"},
		{"chrUn_b", "Un_b"},
		{"chrUn_c", "Un_c"},
		{"chrUn_d", "Un_d"},
		{"chrUn_e", "Un_e"},
	} {
		if !strings.Contains(header, "##contig=<ID="+tc.want+",") {
			t.Errorf("contig %s was not renamed to %s:\n%s", tc.old, tc.want, header)
		}
		if strings.Contains(header, "##contig=<ID="+tc.old+",") {
			t.Errorf("contig %s kept its UCSC spelling:\n%s", tc.old, header)
		}
	}
	// Nothing was lost or duplicated on the way through.
	if got := strings.Count(header, "##contig="); got != 7 {
		t.Errorf("%d contig lines, want 7:\n%s", got, header)
	}
}

func TestStripDropsEveryNamedDefinition(t *testing.T) {
	// --info/--format/--filter are repeatable, not comma-separated.
	out := runVcf(t, "vcf-strip",
		"--info", "IA", "--info", "IB", "--info", "IC", "--info", "ID_",
		"--format", "FA", "--format", "FB", "--format", "FC", "--format", "FD",
		"--filter", "fa", "--filter", "fb", "--filter", "fc", "--filter", "fd",
		"testdata/manycontigs.vcf")

	header := headerOf(out)
	for _, id := range []string{"IA", "IB", "IC", "ID_"} {
		if strings.Contains(header, "##INFO=<ID="+id+",") {
			t.Errorf("##INFO %s survived:\n%s", id, header)
		}
		if strings.Contains(out, id+"=") {
			t.Errorf("INFO key %s survived in the records:\n%s", id, out)
		}
	}
	for _, id := range []string{"FA", "FB", "FC", "FD"} {
		if strings.Contains(header, "##FORMAT=<ID="+id+",") {
			t.Errorf("##FORMAT %s survived:\n%s", id, header)
		}
	}
	for _, id := range []string{"fa", "fb", "fc", "fd"} {
		if strings.Contains(header, "##FILTER=<ID="+id+",") {
			t.Errorf("##FILTER %s survived:\n%s", id, header)
		}
	}

	// The unnamed INFO field and the GT column survive, so this is a projection
	// rather than a wipe.
	if !strings.Contains(header, "##INFO=<ID=KEEPME,") {
		t.Errorf("an INFO field nobody asked to strip was removed:\n%s", header)
	}
	if !strings.Contains(out, "KEEPME=9") {
		t.Errorf("KEEPME was stripped from the records:\n%s", out)
	}
	if !strings.Contains(header, "##FORMAT=<ID=GT,") {
		t.Errorf("GT was removed:\n%s", header)
	}
	// Every remaining sample column is bare GT, with none of the stripped keys.
	for _, line := range strings.Split(out, "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		cols := strings.Split(line, "\t")
		if len(cols) < 10 {
			t.Errorf("record lost its sample columns: %q", line)
			continue
		}
		if cols[8] != "GT" {
			t.Errorf("FORMAT column = %q, want GT: %q", cols[8], line)
		}
	}
}

// headerOf returns just the "##" lines, so a failure message shows the header
// rather than the whole file.
func headerOf(out string) string {
	var b strings.Builder
	for _, line := range strings.Split(out, "\n") {
		if strings.HasPrefix(line, "##") {
			b.WriteString(line + "\n")
		}
	}
	return b.String()
}
