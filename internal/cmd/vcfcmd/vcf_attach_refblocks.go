package vcfcmd

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/compgenlab/cghts/iosource"
	"github.com/compgenlab/cghts/varstore"
	"github.com/compgenlab/cgkit/internal/buildinfo"
	"github.com/spf13/cobra"
)

var (
	attachRebuild bool
	attachReplace []string
	attachFrom    string
	attachBatch   int
	attachTempDir string
)

var vcfVarstoreAttachCmd = &cobra.Command{
	Annotations: map[string]string{"since": "v0.6.11"},
	Use:         "attach-refblocks <store> [fragment ...]",
	Short:       "Merge reference-block fragments into a store's ref_blocks",
	Long: `Merge the per-sample fragments vcf-refblocks wrote into a store's ref_blocks
table, so the store answers at positions off its catalog.

  cgkit vcf-varstore attach-refblocks cohort/ frags/
  cgkit vcf-varstore attach-refblocks cohort/ --from fragments.txt
  cgkit vcf-varstore attach-refblocks cohort/ frags/S0042 --replace S0042

A fragment argument is a fragment directory, or a local directory whose
subdirectories are fragments; --from reads one per line, for when there are
too many for a command line. Every fragment must be finished, one per sample,
for a sample on the store's roster, and built with the store's --min-dp and
depth bands; fragments must also agree with each other on --max-gap.

Blocks are attached only for chromosomes the store's volumes hold, which is the
chromosomes with at least one catalog site. A fragment is genome-wide, and
every volume reads only its own chromosomes from it.

The first attach builds the table from the fragments. After that, say which:

  --rebuild           build it again from these fragments alone; the old
                      table is ignored. The fragments are the source of truth.
  --replace S1,S2     keep the existing table for everyone else, and swap in
                      these samples' fragments -- for a patient whose gVCF was
                      redone. A replaced sample with no fragment is removed.
                      Unchanged fragments need not have been kept.

Either way every file a changed sample appears in is rewritten, because the
table is sorted by position across all samples. What differs is what is read.

Fragments are merged in batches of --batch (default 100), through local
intermediate files, so a ten-thousand-sample cohort never holds ten thousand
files open. The output is a new generation of files and each volume's manifest
is rewritten last: until then the volume reads exactly as before, and a run
that dies leaves the store as it was. On an archive every volume is attached
and the archive manifest rewritten after them; a run that dies partway leaves
volumes that disagree, which the archive refuses to open until it is rerun.`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		ctx := cmd.Context()
		base := strings.TrimSuffix(args[0], "/")

		dirs, err := expandFragments(args[1:], attachFrom)
		if err != nil {
			return err
		}
		if len(dirs) == 0 && len(attachReplace) == 0 {
			return fmt.Errorf("no fragments given")
		}
		frags, err := varstore.LoadFragments(ctx, dirs)
		if err != nil {
			return err
		}
		opts := varstore.AttachOpts{
			Rebuild: attachRebuild,
			Replace: attachReplace,
			Batch:   attachBatch,
			TempDir: attachTempDir,
		}

		volumes := []string{base}
		archive := varstore.IsStore(ctx, base)
		var names []string
		if archive {
			arc, err := varstore.OpenStore(ctx, base)
			if err != nil {
				return err
			}
			volumes = volumes[:0]
			for _, v := range arc.StoreManifest().Volumes {
				names = append(names, v.Name)
				volumes = append(volumes, volumeLocator(base, v.Name))
			}
			arc.Close()
		}

		out := cmd.OutOrStdout()
		fmt.Fprintf(out, "attaching %d fragment(s) to %d volume(s)\n", len(frags), len(volumes))
		for _, v := range volumes {
			res, err := varstore.AttachRefBlocks(ctx, v, frags, opts)
			if err != nil {
				return fmt.Errorf("%s: %w", v, err)
			}
			fmt.Fprintf(out, "  %-32s %12s blocks in %d file(s), %d sample(s)\n",
				v, comma(res.Rows), res.Shards, res.Samples)
		}
		if archive {
			// Every volume now carries ref_blocks; the archive's own manifest
			// still says it does not, and the two must agree before it opens.
			man, err := varstore.BuildStore(ctx, base, names, buildinfo.String(), buildinfo.CommandLine())
			if err != nil {
				return err
			}
			if err := varstore.WriteStoreManifest(base, *man); err != nil {
				return err
			}
		}
		return nil
	},
}

// volumeLocator resolves an archive's volume name the way the archive does: a
// bare name is relative to the store, a locator with a scheme or an absolute
// path is used as given.
func volumeLocator(base, name string) string {
	if strings.Contains(name, "://") || strings.HasPrefix(name, "/") {
		return name
	}
	return base + "/" + name
}

// expandFragments turns arguments and a --from list into fragment directories.
// A local directory without a fragment.json is taken as a parent of fragments,
// which is how ten thousand of them arrive.
func expandFragments(args []string, from string) ([]string, error) {
	var out []string
	add := func(a string) error {
		if iosource.IsRemote(a) || hasFragmentMeta(a) {
			out = append(out, a)
			return nil
		}
		entries, err := os.ReadDir(a)
		if err != nil {
			return fmt.Errorf("%s is neither a fragment nor a directory of them: %w", a, err)
		}
		var found []string
		for _, e := range entries {
			if e.IsDir() && hasFragmentMeta(filepath.Join(a, e.Name())) {
				found = append(found, filepath.Join(a, e.Name()))
			}
		}
		if len(found) == 0 {
			return fmt.Errorf("%s holds no finished fragments", a)
		}
		sort.Strings(found)
		out = append(out, found...)
		return nil
	}
	for _, a := range args {
		if err := add(a); err != nil {
			return nil, err
		}
	}
	if from != "" {
		f, err := os.Open(from)
		if err != nil {
			return nil, err
		}
		defer f.Close()
		sc := bufio.NewScanner(f)
		for sc.Scan() {
			line := strings.TrimSpace(sc.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			out = append(out, line)
		}
		if err := sc.Err(); err != nil {
			return nil, err
		}
	}
	return out, nil
}

func hasFragmentMeta(dir string) bool {
	_, err := os.Stat(filepath.Join(dir, varstore.FragmentMetaFile))
	return err == nil
}

func init() {
	f := vcfVarstoreAttachCmd.Flags()
	f.BoolVar(&attachRebuild, "rebuild", false, "Rebuild ref_blocks from these fragments alone")
	f.StringSliceVar(&attachReplace, "replace", nil, "Samples to replace in the existing ref_blocks")
	f.StringVar(&attachFrom, "from", "", "Read fragment directories from FILE, one per line")
	f.IntVar(&attachBatch, "batch", varstore.DefaultAttachBatch, "Fragments merged at once")
	f.StringVar(&attachTempDir, "tmp", "", "Directory for intermediate merge files")
}
