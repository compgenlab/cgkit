// Package locator decides what a command-line path means and which transports
// this binary can speak.
//
// It exists as a leaf package, rather than as a helper inside internal/cmd,
// for two reasons. The command packages need it and cannot import cmd, since
// root.go imports them. And the blank import below is the only thing that makes
// s3:// work anywhere in cgkit -- putting it in main.go would leave it unlinked
// during `go test ./internal/...`, so a test asserting the "no transport
// registered (registered: s3)" message would pass while the shipped binary
// said something else. Here, that test is the guard on the import.
package locator

import (
	"fmt"

	"github.com/compgenlab/cghts/iosource"

	// Registers the s3:// transport. This is why cgkit links the AWS SDK.
	// iosource.Register panics on a duplicate scheme, but Go runs a package's
	// init once however many times it is imported, so additional blank imports
	// elsewhere are harmless and need no guard.
	_ "github.com/compgenlab/cghts/iosource/s3"
)

// IsRemote reports whether a value names something other than a local file.
//
// It wraps iosource so that cgkit classifies a locator exactly the way the
// library will dispatch it. The rule is not obvious -- a single character
// before "://" is a Windows drive letter, not a scheme -- and a classifier that
// disagreed with the opener by one case would not fail loudly. It would decide
// a URL was a filename, or stat something it should have fetched.
func IsRemote(s string) bool { return iosource.IsRemote(s) }

// Schemes lists the remote transports this binary speaks, for help text.
func Schemes() []string { return append([]string{"http", "https"}, iosource.Schemes()...) }

// CheckLocalOutput rejects a remote locator given as an output destination.
//
// cgkit reads from anywhere and writes locally. Writing to an object store is
// not a small extension of this: there is no writer-side equivalent of a byte
// source, and a multipart upload is a different shape of thing from an
// io.Writer that a BAM or Parquet writer can stream into. Refusing early and by
// name is better than failing deep inside a writer with a mangled path.
//
// Callers must apply this *after* their own stdout check, or "-" gets rejected.
func CheckLocalOutput(flag, value string) error {
	if !IsRemote(value) {
		return nil
	}
	return fmt.Errorf("%s: %q is a remote locator; cgkit writes output to local files only "+
		"(or - for stdout)", flag, value)
}
