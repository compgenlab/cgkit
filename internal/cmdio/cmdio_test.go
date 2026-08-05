package cmdio

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/spf13/cobra"
)

func testCmd() (*cobra.Command, *bytes.Buffer) {
	c := &cobra.Command{Use: "x"}
	var buf bytes.Buffer
	c.SetOut(&buf)
	return c, &buf
}

// Stdout goes through cmd.OutOrStdout(), not os.Stdout. That is the difference
// between a command a test harness can drive and one it cannot, which is why
// several commands had no tests before this package existed.
func TestStdoutGoesThroughTheCommand(t *testing.T) {
	for _, value := range []string{"", "-"} {
		cmd, buf := testCmd()
		out, err := CreateOutput(cmd, "-o/--output", value)
		if err != nil {
			t.Fatalf("CreateOutput(%q): %v", value, err)
		}
		fmt.Fprint(out.W, "hello")
		if err := out.Close(); err != nil {
			t.Errorf("Close on stdout: %v", err)
		}
		if buf.String() != "hello" {
			t.Errorf("value %q wrote %q, want it captured by the command", value, buf.String())
		}
	}
}

func TestFileOutputRoundTrips(t *testing.T) {
	path := filepath.Join(t.TempDir(), "out.txt")
	cmd, buf := testCmd()
	out, err := CreateOutput(cmd, "-o/--output", path)
	if err != nil {
		t.Fatal(err)
	}
	fmt.Fprint(out.W, "hello")
	if err := out.Close(); err != nil {
		t.Fatal(err)
	}
	if buf.Len() != 0 {
		t.Errorf("a file destination also wrote to stdout: %q", buf.String())
	}
	got, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != "hello" {
		t.Errorf("file holds %q, want %q", got, "hello")
	}
}

// The usual pattern is defer Close for the error paths plus an explicit Close at
// the end whose error is returned, so the second call must be harmless.
func TestCloseIsSafeTwiceAndOnStdout(t *testing.T) {
	path := filepath.Join(t.TempDir(), "out.txt")
	cmd, _ := testCmd()
	out, err := CreateOutput(cmd, "-o/--output", path)
	if err != nil {
		t.Fatal(err)
	}
	if err := out.Close(); err != nil {
		t.Fatalf("first Close: %v", err)
	}
	if err := out.Close(); err != nil {
		t.Errorf("second Close: %v, want nil", err)
	}

	stdout, err := CreateOutput(cmd, "-o/--output", "-")
	if err != nil {
		t.Fatal(err)
	}
	if err := stdout.Close(); err != nil {
		t.Errorf("Close on stdout: %v", err)
	}
	// Closing stdout for real would take the process's own stdout with it.
	if stdout.closer != nil {
		t.Error("a stdout destination should own no closer")
	}
}

func TestDiscardRemovesTheFile(t *testing.T) {
	path := filepath.Join(t.TempDir(), "out.txt")
	cmd, _ := testCmd()
	out, err := CreateOutput(cmd, "-o/--output", path)
	if err != nil {
		t.Fatal(err)
	}
	fmt.Fprint(out.W, "partial")
	out.Discard()
	if _, err := os.Stat(path); err == nil {
		t.Error("Discard left the file behind")
	}
}

// Discard on stdout must not panic or try to unlink anything -- the bytes are
// already gone and there is no name to remove.
func TestDiscardOnStdoutIsANoOp(t *testing.T) {
	cmd, buf := testCmd()
	out, err := CreateOutput(cmd, "-o/--output", "-")
	if err != nil {
		t.Fatal(err)
	}
	fmt.Fprint(out.W, "already sent")
	out.Discard()
	if buf.String() != "already sent" {
		t.Errorf("stdout content changed: %q", buf.String())
	}
}

// cgkit reads from anywhere and writes locally, and the refusal names the flag
// it came from rather than failing deep inside a writer with a mangled path.
func TestRemoteOutputIsRefusedByFlagName(t *testing.T) {
	cmd, _ := testCmd()
	_, err := CreateOutput(cmd, "--stats", "s3://bucket/out.txt")
	if err == nil {
		t.Fatal("a remote locator was accepted as an output")
	}
	if !strings.Contains(err.Error(), "--stats") {
		t.Errorf("error %q should name the flag it came from", err)
	}
	if !strings.Contains(err.Error(), "remote locator") {
		t.Errorf("error %q should explain the refusal", err)
	}
}

func TestCreateOutputReportsAnUnwritablePath(t *testing.T) {
	cmd, _ := testCmd()
	dir := t.TempDir()
	// A path whose parent is a regular file cannot be created.
	blocker := filepath.Join(dir, "file")
	if err := os.WriteFile(blocker, nil, 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := CreateOutput(cmd, "-o/--output", filepath.Join(blocker, "out.txt")); err == nil {
		t.Error("expected an error creating a file under a regular file")
	}
}
