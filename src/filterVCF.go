// This scrpt contains functions for applying filters from the interval struct to a vcf file

package main

import (
	"fmt"
	"github.com/icwells/go-tools/iotools"
	"os"
	"strconv"
	"strings"
)

/*func filterNAB(reg Intervals, row []string) bool {
	// Returns true if row passes filtering for coverage in normal file

}*/

func filterCovB(reg Intervals, row []string) bool {
	// Returns true if row passes filtering for coverage in B
	ret := false
	p := reg.Settings
	loc, err := strings.ParseInt(row[1], 10, 64)
	if err == nil {
		_, exists := reg.Regions[row[0]]
		if exists == true {
			_, e := reg.Regions[row[0]][loc]
			if e == true {
				ref = true
				alt = true
				prop = true
				r := reg.Regions[row[0]][loc]
				// Store false if any filter is greater than 0 and is not met
				if p.minb > 0 && r.RefReads <= p.minb {
					ref = false
				}
				if p.maxaltb > 0 && r.AltReads >= p.maxaltb {
					alt = false
				}
				if p.maxprob > 0.0 && (r.AltReads/r.RefReads) >= p.maxprob {
					prop = false
				}
				if ref == true && alt == true && prop == true {
					// Store true if all filters are passed
					ret = true
				}
			}
		}
	}
	return ret
}

func filterVCF(reg Intervals, infile, outfile string) {
	// Applys filters stored in reg to input vcf and writes passing variants to outfile
	fmt.Printf("\n\tFiltering vcf: %s\n", infile)
	f := iotools.OpenFile(infile)
	out := iotools.CreateFile(outfile)
	defer f.Close()
	defer out.Close()
	input := iotools.GetScanner(f)
	for input.Scan() {
		if strings.Contains(line, "#") == true {
			// Write header unchanged
			err := out.WriteString(line)
			_ := iotools.CheckError("Writing header line", err, 0)
		} else {
			var res bool
			s := strings.Split(line, "\t")
			// Apply approriate filters
			if reg.Settings.nab == false {
				res = filterCovB(reg, s)
				//} else {
				//res = filterNAB(reg, s)
			}
			if res == true {
				// Write lines that passed filtering
				err := out.WriteString(line)
				_ := iotools.CheckError("Writing filtered line", err, 0)
			}
		}
	}
}
