package bedcmd

import "testing"

func TestBedToFastaDefault(t *testing.T) {
	bedToFastaIgnoreStrand = false // package-global flag; reset for test ordering
	// Default is stranded: the minus-strand region is reverse-complemented and
	// its coords reported high-to-low. Matches bare `ngsutilsj bed-tofasta`.
	want := ">plusReg|chr1:0-8\nACGTACGT\n" +
		">minusReg|chr1:48-40\nTTTTGGGG\n" +
		">noneReg|chr2:0-4\nTTTT\n" +
		">|chr1:10-14\nGTAA\n"
	if got := runBed(t, "bed-tofasta", "testdata/tofasta.bed", "testdata/ref.fa"); got != want {
		t.Errorf("bed-tofasta default mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestBedToFastaIgnoreStrand(t *testing.T) {
	// --ignore-strand returns plus-strand bases for every region, coords as
	// stored. Matches `ngsutilsj bed-tofasta --ns`.
	want := ">plusReg|chr1:0-8\nACGTACGT\n" +
		">minusReg|chr1:40-48\nCCCCAAAA\n" +
		">noneReg|chr2:0-4\nTTTT\n" +
		">|chr1:10-14\nGTAA\n"
	if got := runBed(t, "bed-tofasta", "--ignore-strand", "testdata/tofasta.bed", "testdata/ref.fa"); got != want {
		t.Errorf("bed-tofasta --ignore-strand mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestBedToFastaWrap(t *testing.T) {
	bedToFastaIgnoreStrand = false // package-global flag; reset for test ordering
	// Default (stranded) with wrapping at 4 bases.
	want := ">plusReg|chr1:0-8\nACGT\nACGT\n" +
		">minusReg|chr1:48-40\nTTTT\nGGGG\n" +
		">noneReg|chr2:0-4\nTTTT\n" +
		">|chr1:10-14\nGTAA\n"
	if got := runBed(t, "bed-tofasta", "--wrap", "4", "testdata/tofasta.bed", "testdata/ref.fa"); got != want {
		t.Errorf("bed-tofasta --wrap mismatch.\n got: %q\nwant: %q", got, want)
	}
}
