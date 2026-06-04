package cmd

import "testing"

// TestInit forces the cmd package (and all subcommand packages it imports) to
// load under `go test ./...`. Cobra/pflag panic at init() on conditions like
// duplicate flag names, so loading the package is enough to surface those —
// without this test, no test file in the binary's package tree would trigger
// the init() chain and these panics would only show up in built binaries.
func TestInit(t *testing.T) {
	if rootCmd == nil {
		t.Fatal("rootCmd is nil")
	}
}
