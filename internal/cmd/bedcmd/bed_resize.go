package bedcmd

import (
	"fmt"
	"io"

	"github.com/compgenlab/hts/bed"
	"github.com/spf13/cobra"
)

var (
	bedResizeOutput string
	bedResizeLen5   int
	bedResizeLen3   int
	bedResizeMax    int
)

var bedResizeCmd = &cobra.Command{
	GroupID:     "bedcmd",
	Annotations: map[string]string{"since": "v0.3.1"},
	Use:         "bed-resize <input.bed>",
	Short:       "Resize BED regions (extend or shrink)",
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) == 0 {
			cmd.Help()
			return nil
		}
		if bedResizeLen5 == 0 && bedResizeLen3 == 0 {
			return fmt.Errorf("you must set either -5 or -3 (or both)")
		}

		reader, err := openBedInput(cmd, args[0])
		if err != nil {
			return err
		}
		defer reader.Close()

		// Resized regions retain name/score/strand/extras and use float scores,
		// matching ngsutilsj bed-resize.
		writer, err := openBedOutput(cmd, bedResizeOutput, bed.NewBedWriterOpts())
		if err != nil {
			return err
		}

		for {
			rec, err := reader.NextRecord()
			if err != nil {
				if err == io.EOF {
					break
				}
				return err
			}

			// Regions already at or above the max length are passed through
			// unchanged.
			if bedResizeMax <= 0 || rec.Length() < bedResizeMax {
				rec = rec.Extend5(bedResizeLen5).Extend3(bedResizeLen3)
			}

			if err := writer.WriteRecord(rec); err != nil {
				return err
			}
		}
		return writer.Close()
	},
}

func init() {
	bedResizeCmd.Flags().IntVarP(&bedResizeLen5, "extend5", "5", 0, "Extend a region in the 5' direction (strand-specific, negative to shrink)")
	bedResizeCmd.Flags().IntVarP(&bedResizeLen3, "extend3", "3", 0, "Extend a region in the 3' direction (strand-specific, negative to shrink)")
	bedResizeCmd.Flags().IntVar(&bedResizeMax, "max", -1, "Maximum length to expand (regions at or above this length are not expanded)")
	bedResizeCmd.Flags().StringVarP(&bedResizeOutput, "output", "o", "-", "Output filename (gzip-compressed if it ends in .gz; - for stdout)")
}
