// Package cmdio opens command inputs and outputs the same way everywhere.
//
// It exists because the same twelve lines -- is it stdout, is it a remote
// locator, create it, remember to close it -- were written out by hand in five
// command packages, and the copies had drifted. Three defaulted to
// cmd.OutOrStdout() and two to os.Stdout, which is not a style difference: a
// command writing to os.Stdout cannot be driven by a cobra test harness at all,
// and that is most of why several commands had no tests. Three returned the
// close error and two discarded it, so a full disk produced a truncated file and
// a successful exit.
//
// Like internal/locator, this is a leaf package: the command packages need it
// and cannot import internal/cmd, since root.go imports them.
package cmdio

import (
	"fmt"
	"io"
	"os"

	"github.com/compgenlab/cgkit/internal/locator"
	"github.com/spf13/cobra"
)

// Output is a command's destination, and whether closing it is the caller's
// problem.
type Output struct {
	// W is where to write. Never nil.
	W io.Writer

	// closer is nil for stdout, which this package does not own and must not
	// close -- doing so would take the process's stdout down with it.
	closer io.Closer
	name   string
}

// Close finishes the output and reports any error.
//
// It is safe on a stdout destination and safe to call twice, so the usual
// pattern works: defer it for the error paths, and call it explicitly at the end
// where the error can actually be returned. Discarding it is how a failed flush
// became a successful exit in three commands.
func (o *Output) Close() error {
	if o == nil || o.closer == nil {
		return nil
	}
	c := o.closer
	o.closer = nil
	if err := c.Close(); err != nil {
		return fmt.Errorf("closing %s: %w", o.name, err)
	}
	return nil
}

// Discard closes the output and removes it, for a run that failed partway.
//
// A half-written file is worse than none: it looks like a result, and for a
// compressed format it is a stream with no terminator, which is detectably
// broken only if something checks. Stdout is left alone -- there is nothing to
// remove and the bytes are already gone.
func (o *Output) Discard() {
	if o == nil || o.name == "" {
		return
	}
	o.Close()
	os.Remove(o.name)
}

// CreateOutput opens a command's output. An empty value or "-" means stdout.
//
// flag names the flag being resolved, so a rejected remote locator can say which
// one it came from.
func CreateOutput(cmd *cobra.Command, flag, value string) (*Output, error) {
	if value == "" || value == "-" {
		return &Output{W: cmd.OutOrStdout()}, nil
	}
	// After the stdout branch, never before: "-" is not a locator. Harmless
	// either way now that CheckLocalOutput accepts stdout itself, but the order
	// still says what is meant.
	if err := locator.CheckLocalOutput(flag, value); err != nil {
		return nil, err
	}
	f, err := os.Create(value)
	if err != nil {
		return nil, err
	}
	return &Output{W: f, closer: f, name: value}, nil
}
