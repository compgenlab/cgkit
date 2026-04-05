package fastqcmd

import (
	"bytes"
	"os"
	"testing"

	"github.com/spf13/cobra"
)

func TestFastqTag(t *testing.T) {
	// Write a test FASTQ file
	tmp, err := os.CreateTemp("", "test_tag_*.fq")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmp.Name())
	tmp.WriteString("@read1\nACGT\n+\nIIII\n@read2 existing\nTTTT\n+\nFFFF\n")
	tmp.Close()

	root := &cobra.Command{}
	root.AddGroup(&cobra.Group{ID: "fastqcmd", Title: "FASTQ"})
	root.AddCommand(fastqTagCmd)

	buf := new(bytes.Buffer)
	root.SetOut(buf)
	root.SetArgs([]string{"fastq-tag", "BC:Z:ACGT", tmp.Name()})

	if err := root.Execute(); err != nil {
		t.Fatal(err)
	}

	got := buf.String()
	expected := "@read1 BC:Z:ACGT\nACGT\n+\nIIII\n@read2 existing\tBC:Z:ACGT\nTTTT\n+\nFFFF\n"
	if got != expected {
		t.Errorf("got:\n%s\nexpected:\n%s", got, expected)
	}
}
