package main

import (
	"reflect"
	"testing"
)

// os.Args is provenance: buildinfo.CommandLine joins it into the ##...Command
// header of every VCF, the @PG CL: of every BAM and cgkit.command in every
// store. The previous code re-sliced os.Args[1:], which drops the *program
// name* rather than the flag -- so every record written under --profile was
// stamped "--profile=cpu.prof vcf-stats in.vcf", with no "cgkit" in front.
func TestTakeProfileFlag(t *testing.T) {
	cases := []struct {
		name     string
		in       []string
		wantArgs []string
		wantPath string
	}{
		{
			name:     "absent",
			in:       []string{"cgkit", "vcf-stats", "in.vcf"},
			wantArgs: []string{"cgkit", "vcf-stats", "in.vcf"},
		},
		{
			name:     "before the subcommand",
			in:       []string{"cgkit", "--profile=cpu.prof", "vcf-stats", "in.vcf"},
			wantArgs: []string{"cgkit", "vcf-stats", "in.vcf"},
			wantPath: "cpu.prof",
		},
		{
			// The order anyone who has used a CLI would reach for. It used to
			// die with "unknown flag" because only os.Args[1] was examined.
			name:     "after the subcommand",
			in:       []string{"cgkit", "vcf-stats", "--profile=cpu.prof", "in.vcf"},
			wantArgs: []string{"cgkit", "vcf-stats", "in.vcf"},
			wantPath: "cpu.prof",
		},
		{
			name:     "at the end",
			in:       []string{"cgkit", "vcf-stats", "in.vcf", "--profile=cpu.prof"},
			wantArgs: []string{"cgkit", "vcf-stats", "in.vcf"},
			wantPath: "cpu.prof",
		},
		{
			name:     "empty path is still consumed",
			in:       []string{"cgkit", "--profile=", "vcf-stats"},
			wantArgs: []string{"cgkit", "vcf-stats"},
			wantPath: "",
		},
		{
			// Everything past a bare "--" is positional by convention.
			name:     "after a bare double dash",
			in:       []string{"cgkit", "vcf-stats", "--", "--profile=cpu.prof"},
			wantArgs: []string{"cgkit", "vcf-stats", "--", "--profile=cpu.prof"},
		},
		{
			// A filename that merely starts the same way is not the flag.
			name:     "not a prefix match on a positional",
			in:       []string{"cgkit", "vcf-stats", "--profiles=x"},
			wantArgs: []string{"cgkit", "vcf-stats", "--profiles=x"},
		},
		{
			// argv[0] is never a flag, however it is spelled.
			name:     "argv0 is never consumed",
			in:       []string{"--profile=weird-name"},
			wantArgs: []string{"--profile=weird-name"},
		},
		{
			name:     "only the first is consumed",
			in:       []string{"cgkit", "--profile=a.prof", "vcf-stats", "--profile=b.prof"},
			wantArgs: []string{"cgkit", "vcf-stats", "--profile=b.prof"},
			wantPath: "a.prof",
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			args, path := takeProfileFlag(c.in)
			if !reflect.DeepEqual(args, c.wantArgs) {
				t.Errorf("args = %q, want %q", args, c.wantArgs)
			}
			if path != c.wantPath {
				t.Errorf("path = %q, want %q", path, c.wantPath)
			}
			// The program name has to survive, or every provenance record
			// written by the run is missing the tool that wrote it.
			if len(args) == 0 || args[0] != c.in[0] {
				t.Errorf("args[0] = %q, want the program name %q", args, c.in[0])
			}
		})
	}
}
