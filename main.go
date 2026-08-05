package main

import (
	"fmt"
	"os"
	"runtime/pprof"
	"strings"

	"github.com/compgenlab/cgkit/internal/cmd"
)

// profileFlag is intercepted here rather than declared on the root command,
// because profiling has to start before cobra dispatches and stop after it
// returns.
const profileFlag = "--profile="

func main() {
	os.Exit(run())
}

// run is separate from main so that its defers execute. os.Exit does not run
// them, and the whole point of the profile teardown is that it happens even
// when the command fails -- a failing run used to leave a zero-byte, unparseable
// profile because the deferred StopCPUProfile never fired.
func run() int {
	args, path := takeProfileFlag(os.Args)
	// os.Args is provenance: buildinfo.CommandLine joins it into the ##...Command
	// header of every VCF, the @PG CL: of every BAM and cgkit.command in every
	// store. It is rewritten rather than left alone so those records describe a
	// reproducible invocation instead of one carrying a profiling flag -- and the
	// program name is kept, which the previous re-slicing dropped, stamping
	// "--profile=cpu.prof vcf-stats in.vcf" with no "cgkit" in front of it.
	os.Args = args

	if path != "" {
		f, err := os.Create(path)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error: --profile: %v\n", err)
			return 1
		}
		defer func() {
			pprof.StopCPUProfile()
			if err := f.Close(); err != nil {
				fmt.Fprintf(os.Stderr, "Error: --profile: %v\n", err)
			}
		}()
		if err := pprof.StartCPUProfile(f); err != nil {
			fmt.Fprintf(os.Stderr, "Error: --profile: %v\n", err)
			f.Close()
			return 1
		}
	}

	if err := cmd.Run(); err != nil {
		// Cobra has already printed it.
		return 1
	}
	return 0
}

// takeProfileFlag removes a --profile=FILE argument from args and returns the
// remainder along with the filename.
//
// It is scanned for anywhere rather than only at args[1]: the positional form
// meant "cgkit vcf-stats --profile=cpu.prof in.vcf" -- the order anyone who has
// used a CLI would reach for -- died with "unknown flag" instead of profiling.
// Everything after a bare "--" is a positional argument by convention and is
// left alone.
func takeProfileFlag(args []string) (rest []string, path string) {
	rest = make([]string, 0, len(args))
	for i, a := range args {
		if i > 0 && a == "--" {
			rest = append(rest, args[i:]...)
			break
		}
		if v, ok := strings.CutPrefix(a, profileFlag); ok && i > 0 && path == "" {
			path = v
			continue
		}
		rest = append(rest, a)
	}
	return rest, path
}
