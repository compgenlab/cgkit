package main

import (
	"log"
	"os"
	"runtime/pprof"

	"github.com/compgen-io/cgltk/internal/cmd"
)

func main() {
	if len(os.Args) > 1 {
		if os.Args[1][0:10] == "--profile=" {
			pfile := os.Args[1][10:len(os.Args[1])]
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
