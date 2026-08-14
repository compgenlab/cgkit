package vcfcmd

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/compgenlab/cghts/varstore"
	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/spf13/cobra"
)

// A varstore: several volumes, disjoint by chromosome, read as one.
//
// WHY THIS IS A SEPARATE STEP rather than something a conversion does. A
// whole-genome callset arrives as one VCF per chromosome and is converted one
// volume at a time -- often in parallel, often on different machines, often
// days apart. Combining is the statement that those volumes belong together,
// and it can only be made once they all exist, because making it means checking
// that they agree.
//
// It also means a varstore can be rebuilt, or one volume replaced, without
// reconverting anything.
//
// Note that vcf-tovarstore already produces a varstore -- a complete one of a
// single volume. This does not promote anything to a different kind of object;
// it composes volumes into a larger archive of the same kind.

var vcfVarstoreCmd = &cobra.Command{
	GroupID:     "vcfcmd",
	Annotations: map[string]string{"since": "v0.6.9"},
	Use:         "vcf-varstore",
	Short:       "Combine volumes into one varstore, queried as one",
	Long: `A varstore is several volumes, disjoint by chromosome, read as one.

Whole-genome data ships one VCF per chromosome, so converting it gives one
volume per chromosome -- and twenty-five volumes is the wrong unit to manage,
register or reason about. A person has one genome; a varstore is one object that
means it.

It is also what makes small shards affordable. A shard index lives in a manifest
read on every open, so its size is a budget: chr22 split every 500 sites carries
about 100 KB of index, while a whole genome split the same way carries 2.6 MB. A
single undivided archive forces a choice between coarse shards and an unreadable
manifest, and coarse shards are what make a locus query slow. Split by
chromosome, each index covers one chromosome and a query opens only the
chromosomes it touches.

The words, since two of them used to be one word:

  varstore   the whole archive -- what you hold, name and query
  volume     one chromosome's worth of it, a complete varstore on its own
  shard      a coordinate range within a volume`,
}

var vcfVarstoreCreateCmd = &cobra.Command{
	// No GroupID: groups belong to the root's listing, and a subcommand of a
	// subcommand is not in it.
	Annotations: map[string]string{"since": "v0.6.9"},
	Use:         "create <store-dir> [volume ...]",
	Short:       "Write a store manifest over volumes that already exist",
	Long: `Combine existing volumes into one varstore.

  cgkit vcf-varstore create cohort/ chr1 chr2 chr3
  cgkit vcf-varstore create s3://bucket/cohort chr1 chr2
  cgkit vcf-varstore create cohort/

A volume named without a scheme or a leading slash is relative to the store
directory, which keeps a varstore one thing to copy, move or delete. Named with
a locator, it is used as given, so volumes may live apart when they must.

Given no volumes, every subdirectory of the store that is a readable volume is
taken, in name order.

THE VOLUMES MUST AGREE, and this is where that is established rather than
rediscovered on every query:

  the same samples, in the same order    genotype columns are positional, and a
                                         mismatch attributes calls to the wrong
                                         person
  the same --min-dp                      otherwise they do not mean the same
                                         thing by "callable"
  the same --depth-bands                 otherwise a recorded depth describes a
                                         different span in each
  no shared chromosome                   otherwise a locus has two answers and
                                         no rule for choosing

A varstore that cannot be built is a better outcome than one that answers with a
different population depending on which chromosome you ask about.

The manifest is written LAST and named varstore.json.gz. Its presence is what
makes the directory an archive of several volumes, and cgkit and anything else
reading through cghts will then open it as one.`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		base := strings.TrimSuffix(args[0], "/")
		volumes := args[1:]

		if len(volumes) == 0 {
			found, err := discoverVolumes(cmd, base)
			if err != nil {
				return err
			}
			if len(found) == 0 {
				return fmt.Errorf(
					"%s holds no volumes; name them, or convert them into it first", base)
			}
			volumes = found
		}

		man, err := varstore.BuildStore(cmd.Context(), base, volumes,
			buildinfo.String(), buildinfo.CommandLine())
		if err != nil {
			return err
		}
		if err := varstore.WriteStoreManifest(base, *man); err != nil {
			return err
		}

		out := cmd.OutOrStdout()
		fmt.Fprintf(out, "\nwrote %s: %d volumes, %d samples\n\n",
			varstore.StoreManifestPath(base), len(man.Volumes), len(man.Samples))
		for _, m := range man.Volumes {
			fmt.Fprintf(out, "  %-24s %-28s %s sites\n",
				m.Name, strings.Join(m.Chroms, " "), comma(m.Sites))
		}
		fmt.Fprintln(out)
		return nil
	},
}

// discoverVolumes lists the subdirectories of a store that are readable volumes.
//
// Only local paths: listing is not something every locator scheme offers, and a
// caller on object storage knows its own volume names. Naming them is also the
// safer habit -- a discovered archive silently gains whatever was left in the
// directory.
func discoverVolumes(cmd *cobra.Command, base string) ([]string, error) {
	entries, err := osReadDirNames(base)
	if err != nil {
		return nil, fmt.Errorf("listing %s: %w (name the volumes instead)", base, err)
	}
	var out []string
	for _, name := range entries {
		if _, err := varstore.ReadVolumeManifest(joinLocal(base, name)); err != nil {
			continue
		}
		out = append(out, name)
	}
	sort.Strings(out)
	return out, nil
}

func comma(n int64) string {
	s := fmt.Sprint(n)
	if len(s) <= 3 {
		return s
	}
	var b strings.Builder
	for i, r := range s {
		if i > 0 && (len(s)-i)%3 == 0 {
			b.WriteByte(',')
		}
		b.WriteRune(r)
	}
	return b.String()
}

func osReadDirNames(p string) ([]string, error) {
	es, err := os.ReadDir(p)
	if err != nil {
		return nil, err
	}
	var out []string
	for _, e := range es {
		if e.IsDir() {
			out = append(out, e.Name())
		}
	}
	return out, nil
}

func joinLocal(a, b string) string { return filepath.Join(a, b) }
