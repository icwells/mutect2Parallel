// This script will filter vcf files based on output from covB.sh/nab.sh

package main

import (
	"fmt"
	"github.com/icwells/go-tools/iotools"
	"gopkg.in/alecthomas/kingpin.v2"
	"os"
	"strings"
	"time"
)

var (
	app      = kingpin.New("heterAnalyzer", "Applies custom filters to vcf files.")
	infile   = kingpin.Flag("vcf", "Path to input vcf.").Required().Short('v').String()
	interval = kingpin.Flag("interval", "Path to covB.sh/nab.sh output.").Required().Short('i').String()
	outfile  = kingpin.Flag("outfile", "Path to output vcf.").Required().Short('o').String()

	covb     = kingpin.Command("covb", "Filter input based on coverage in vcf from same tumor (B).")
	minb     = covb.Flag("min_covB", "Minimum coverage in B.").Default("15").Int()
	maxaltb  = covb.Flag("max_altB", "Maximum number of alternative alleles in B.").Default("0").Int()
	maxpropb = covb.Flag("max_prop_altB", "Maximum proportion of alternative reads in B.").Default("0.0").Float()

	nab      = kingpin.Command("nab", "Filter based on coverage in normal tissue vcf.")
	minn     = nab.Flag("min_covN", "Minimum coverage in normal vcf.").Default("5").Int()
	maxaltn  = nab.Flag("max_reads_altN", "Maximum number of alternative alleles in normal vcf.").Default("15").Int()
	maxpropn = nab.Flag("max_freq_altN", "Maximum proportion of alternative reads in normal vcf.").Default("0.3").Float()
)

func checkInput(infile string) string {
	// Check for existance of infile with and without ".gz"; exits if it does not exist
	if iotools.Exists(infile) == false {
		if strings.Contains(infile, ".gz") == true {
			if strings.Count(infile, ".") > 1 {
				idx := strings.LastIndex(infile, ".")
				infile = infile[:idx]
			}
		} else {
			infile = infile + ".gz"
		}
	}
	// Repeat with edited infile
	if iotools.Exists(infile) == false {
		fmt.Printf("\n[Error] Cannot find %s. Exiting.\n\n", infile)
		os.Exit(1)
	}
	return infile
}

func main() {
	start := time.Now()
	var reg Intervals
	switch kingpin.Parse() {
	case covb.FullCommand():
		reg.setIntervals(false, *minb, *maxaltb, *minn, *maxaltn, *maxpropb, *maxpropn)
	case nab.FullCommand():
		reg.setIntervals(true, *minb, *maxaltb, *minn, *maxaltn, *maxpropb, *maxpropn)
	}
	i := checkInput(*interval)
	vcf := checkInput(*infile)
	reg.loadIntervals(i)
	filterVCF(reg, vcf, *outfile)
	fmt.Printf("\tFinished. Runtime: %s\n\n", time.Since(start))
}
