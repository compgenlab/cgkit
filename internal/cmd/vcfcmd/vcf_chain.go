package vcfcmd

import "github.com/spf13/pflag"

// chainArg is one chain-flag occurrence: the flag name (kind) and its argument
// (value is "true" for boolean flags).
type chainArg struct {
	kind  string
	value string
}

// chainValue is a pflag.Value that appends to a shared ordered slice when set.
// Because pflag calls Set in command-line order, the slice captures the exact
// order (and interleaving) of the chain flags. Used by vcf-annotate and
// vcf-filter to apply their annotators/filters in command-line order.
type chainValue struct {
	kind   string
	isBool bool
	target *[]chainArg
}

func (c *chainValue) String() string { return "" }

func (c *chainValue) Set(v string) error {
	*c.target = append(*c.target, chainArg{kind: c.kind, value: v})
	return nil
}

func (c *chainValue) Type() string {
	if c.isBool {
		return "bool"
	}
	return "string"
}

// IsBoolFlag makes pflag treat the flag as a no-argument boolean when isBool.
func (c *chainValue) IsBoolFlag() bool { return c.isBool }

// registerChainBool registers a boolean chain flag. NoOptDefVal is required so
// pflag does not consume the next argument as the flag's value.
func registerChainBool(f *pflag.FlagSet, target *[]chainArg, name, usage string) {
	f.Var(&chainValue{kind: name, isBool: true, target: target}, name, usage)
	f.Lookup(name).NoOptDefVal = "true"
}

// registerChainVal registers a value-taking chain flag.
func registerChainVal(f *pflag.FlagSet, target *[]chainArg, name, usage string) {
	f.Var(&chainValue{kind: name, target: target}, name, usage)
}
