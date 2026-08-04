// Package cmdtest holds helpers shared by the command packages' tests.
//
// It is a normal package rather than a _test file because Go cannot share test
// helpers across packages any other way, and nothing outside tests imports it.
package cmdtest

import (
	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// ResetFlags restores every flag on every subcommand of root to its default.
//
// Each test builds a fresh root, but the commands under it are not fresh:
// InitCmd registers package-level cobra.Command values whose flags are bound to
// package-level variables. So a flag set by one test stays set for every test
// that follows, and the failure is order-dependent by construction and never
// loud -- it silently changes some later test's output. bed-tofasta inheriting
// a --wrap from an earlier case is what that looked like in practice, and it
// broke 16 of 25 `go test -shuffle` seeds.
//
// Callers with flags whose values do not round-trip through Set -- custom
// chain-valued flags, mostly -- still have to clear those by hand afterwards.
func ResetFlags(root *cobra.Command) {
	for _, c := range root.Commands() {
		c.Flags().VisitAll(func(f *pflag.Flag) {
			// Only flags something actually set. Changed persists across runs
			// for the same reason the values do -- the flag objects are shared --
			// so it is exactly the record of what needs undoing.
			//
			// Touching the rest is not merely wasted work, it is unsafe. A
			// custom flag type that accumulates has a DefValue rendering as
			// "[]", and its Set accepts that literal happily: sam-filter's --tag
			// then fails with `invalid tag filter "[]"` in a test that never
			// passed --tag at all.
			if !f.Changed {
				return
			}
			// Slice flags APPEND on Set, so Set(DefValue) would add the literal
			// "[]" as an element rather than clearing them. Replace empties them
			// properly, and does so for every slice flag rather than a
			// hand-maintained list that a new flag gets forgotten from.
			if sv, ok := f.Value.(pflag.SliceValue); ok {
				_ = sv.Replace(nil)
				f.Changed = false
				return
			}
			_ = f.Value.Set(f.DefValue)
			f.Changed = false
		})
	}
}
