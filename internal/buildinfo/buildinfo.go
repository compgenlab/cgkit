// Package buildinfo holds the cgio version (set via -ldflags at build time) and
// provenance helpers shared across commands. Keeping it in a leaf package lets
// every command group reference the version without importing the cmd package.
package buildinfo

import (
	"os"
	"strings"
	"time"
)

// Set via -ldflags at build time.
var (
	Version = "dev"
	GitHash = ""
)

// Now returns the current time. It is a package variable so tests can inject a
// fixed clock for deterministic provenance (e.g. VCF ##fileDate) output.
var Now = time.Now

// String returns the version, with the git hash appended when available.
func String() string {
	if GitHash != "" {
		return Version + " (" + GitHash + ")"
	}
	return Version
}

// CommandLine returns the full command-line invocation, for provenance records.
func CommandLine() string {
	return strings.Join(os.Args, " ")
}

// Date returns the current date as YYYYMMDD, the VCF ##fileDate format.
func Date() string {
	return Now().Format("20060102")
}
