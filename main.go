package main

import (
	"log"
	"os"
	"runtime/pprof"
	"strings"

	"github.com/compgenlab/cgio/internal/cmd"
)

func main() {
	if len(os.Args) > 1 {
		if strings.HasPrefix(os.Args[1], "--profile=") {
			pfile := os.Args[1][10:]
			f, err := os.Create(pfile)
			if err != nil {
				log.Fatal(err)
			}
			defer f.Close()

			if err := pprof.StartCPUProfile(f); err != nil {
				log.Fatal(err)
			}
			defer pprof.StopCPUProfile()

			os.Args = os.Args[1:len(os.Args)]
		}
	}

	cmd.Execute()
}
