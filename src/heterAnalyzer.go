// This script will filter vcf files based on output from covB.sh/nab.sh

package main

import (
	"fmt"
	"gopkg.in/alecthomas/kingpin.v2"
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
	maxn     = nab.Flag("max_covN", "Maximum coverage in normal vcf.").Default("5").Int()
	maxaltn  = nab.Flag("min_reads_altN", "Maximum number of alternative alleles in normal vcf.").Default("15").Int()
	maxpropn = nab.Flag("min_freq_altN", "Maximum proportion of alternative reads in normal vcf.").Default("0.3").Float()
)

func main() {
	start := time.Now()
	var reg Intervals
	switch kingpin.Parse() {
	case covb.FullCommand():
		reg.setIntervals(false, *minb, *maxaltb, *maxn, *maxaltn, *maxpropb, *maxpropn)
	case nab.FullCommand():
		reg.setIntervals(true, *minb, *maxaltb, *maxn, *maxaltn, *maxpropb, *maxpropn)
	}
	reg.loadIntervals(*interval)
	filterVCF(reg, *infile, *outfile)
	fmt.Printf("\n\tFinished. Runtime: %s\n\n", time.Since(start))
}
