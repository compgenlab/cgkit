package ontcmd

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/compgen-io/cgltk/htsio"
	"github.com/spf13/cobra"
)

// umiClusterRecord represents one line from the ont-umi-cluster --umi-counts output.
type umiClusterRecord struct {
	chrom          string
	start          int // 0-based
	end            int // 0-based, exclusive
	mi             string
	strand         string
	representative string
	umis           []string
}

// matchesPosition returns true if a read at (start, end) on the given strand
// overlaps this cluster record within the gap tolerance.
func (r *umiClusterRecord) matchesPosition(rname, strand string, readStart, readEnd, gap int, matchOneEnd, noStrand bool) bool {
	if rname != r.chrom {
		return false
	}
	if !noStrand && strand != r.strand {
		return false
	}

	fivePrime := abs(readStart-r.start) <= gap
	threePrime := abs(readEnd-r.end) <= gap

	if matchOneEnd {
		return fivePrime || threePrime
	}
	return fivePrime && threePrime
}

// matchesUMI returns true if the query UMI is within editDist of any UMI in
// this cluster.
func (r *umiClusterRecord) matchesUMI(queryUMI string, maxDist int) bool {
	norm := normalizeUMISeparator(queryUMI)
	var buf levBuf
	for _, umi := range r.umis {
		// Bounded: we only care whether the distance is <= maxDist, so
		// the Ukkonen cutoff lets us bail out early on non-matches.
		if levDist(norm, normalizeUMISeparator(umi), &buf, maxDist) <= maxDist {
			return true
		}
	}
	return false
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// parseCounts reads the UMI cluster counts file (plain text or gzipped).
// Returns records sorted by (chrom, start, end) as written by ont-umi-cluster.
func parseCounts(filename string) ([]umiClusterRecord, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var reader io.Reader = f
	if strings.HasSuffix(filename, ".gz") || strings.HasSuffix(filename, ".bgz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return nil, fmt.Errorf("decompressing %s: %w", filename, err)
		}
		defer gz.Close()
		reader = gz
	}

	var records []umiClusterRecord
	scanner := bufio.NewScanner(reader)
	lineNum := 0
	for scanner.Scan() {
		lineNum++
		line := scanner.Text()
		if line == "" {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 10 {
			return nil, fmt.Errorf("line %d: expected 10 fields, got %d", lineNum, len(fields))
		}

		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("line %d: bad start: %w", lineNum, err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return nil, fmt.Errorf("line %d: bad end: %w", lineNum, err)
		}

		umis := strings.Split(fields[9], ",")

		records = append(records, umiClusterRecord{
			chrom:          fields[0],
			start:          start,
			end:            end,
			mi:             fields[3],
			strand:         fields[4],
			representative: fields[5],
			umis:           umis,
		})
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return records, nil
}

var ontUmiLookupCmd = &cobra.Command{
	GroupID: "ontcmd",
	Use:     "ont-umi-lookup <counts.bed.gz> <query.bam>",
	Short:   "Check if reads in a BAM file match existing UMI clusters",
	Long: `For each read in the query BAM file with a UMI tag, check the
ont-umi-cluster counts file to see if it belongs to an existing UMI cluster
based on positional overlap and UMI edit distance.

The counts file must be the output of ont-umi-cluster --umi-counts.
Both files must be sorted by coordinate.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) < 2 {
			cmd.Help()
			return nil
		}

		countsFile := args[0]
		bamFile := args[1]

		// Load counts records.
		records, err := parseCounts(countsFile)
		if err != nil {
			return fmt.Errorf("reading counts: %w", err)
		}
		fmt.Fprintf(os.Stderr, "Loaded %d UMI cluster records\n", len(records))

		reader, err := htsio.NewSamReader(bamFile)
		if err != nil {
			return err
		}
		defer reader.Close()

		var out io.Writer = os.Stdout
		if umiLookupOutput != "" && umiLookupOutput != "-" {
			f, err := os.Create(umiLookupOutput)
			if err != nil {
				return err
			}
			defer f.Close()
			out = f
		}

		// Write header.
		fmt.Fprintln(out, strings.Join([]string{
			"read_name", "chrom", "start", "end", "strand",
			"umi", "match", "mi", "representative_umi",
		}, "\t"))

		// Sweep through BAM and counts in coordinate order.
		// recordIdx tracks the start of the window of candidate records.
		recordIdx := 0
		matched := 0
		total := 0

		for {
			rec, err := reader.Next()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if rec.IsUnmapped() || rec.Cigar == "*" {
				continue
			}

			umi := ""
			if tag, ok := rec.Tags[umiLookupTag]; ok {
				umi = tag.Value
			}
			if umi == "" {
				continue
			}

			readStart := rec.Pos - 1
			readEnd := readStart + htsio.CigarRefLen(rec.Cigar)
			strand := "+"
			if rec.IsReverse() {
				strand = "-"
			}

			// Advance recordIdx past records that can no longer match
			// any current or future read (sorted input assumption).
			for recordIdx < len(records) {
				r := &records[recordIdx]
				if r.chrom < rec.RefName {
					recordIdx++
					continue
				}
				if r.chrom == rec.RefName && r.end+umiLookupOverlap < readStart {
					recordIdx++
					continue
				}
				break
			}

			// Scan candidate records for matches.
			total++
			bestMI := ""
			bestRep := ""
			for i := recordIdx; i < len(records); i++ {
				r := &records[i]
				// Past this read's possible range — stop scanning.
				if r.chrom > rec.RefName {
					break
				}
				if r.chrom == rec.RefName && r.start-umiLookupOverlap > readStart {
					break
				}

				if !r.matchesPosition(rec.RefName, strand, readStart, readEnd, umiLookupOverlap, umiLookupMatchOneEnd, umiLookupNoStrand) {
					continue
				}
				if r.matchesUMI(umi, umiLookupEditDist) {
					bestMI = r.mi
					bestRep = r.representative
					break
				}
			}

			matchStr := "no"
			if bestMI != "" {
				matchStr = "yes"
				matched++
			}
			fmt.Fprintf(out, "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n",
				rec.ReadName, rec.RefName, readStart, readEnd, strand,
				umi, matchStr, bestMI, bestRep)
		}

		fmt.Fprintf(os.Stderr, "Total reads with UMI: %d, matched: %d\n", total, matched)
		return nil
	},
}

var (
	umiLookupOutput      string
	umiLookupTag         string
	umiLookupOverlap     int
	umiLookupEditDist    int
	umiLookupMatchOneEnd bool
	umiLookupNoStrand    bool
)

func init() {
	ontUmiLookupCmd.Flags().StringVarP(&umiLookupOutput, "output", "o", "-", "Output file (default: stdout)")
	ontUmiLookupCmd.Flags().StringVar(&umiLookupTag, "umi-tag", "RX", "SAM tag containing UMI sequence")
	ontUmiLookupCmd.Flags().IntVar(&umiLookupOverlap, "overlap", 50, "Maximum gap (bp) between read ends to consider a positional match")
	ontUmiLookupCmd.Flags().IntVar(&umiLookupEditDist, "umi-edit-distance", 3, "Maximum Levenshtein edit distance to match a UMI")
	ontUmiLookupCmd.Flags().BoolVar(&umiLookupMatchOneEnd, "match-one-end", false, "Match if EITHER 5' or 3' ends overlap (default: require BOTH)")
	ontUmiLookupCmd.Flags().BoolVar(&umiLookupNoStrand, "no-strand", false, "Ignore strand when matching positions")
}
