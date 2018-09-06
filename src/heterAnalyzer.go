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
	infile   = kingpin.Flag("vcf", "Path to input vcf.").Required().Short('v').Default("nil").String()
	interval = kingpin.Flag("interval", "Path to covB.sh/nab.sh output.").Required().Short('i').Default("nil").String()
	outfile  = kingpin.Flag("outfile", "Path to output vcf.").Required().Short('o').Default("nil").String()

	covb     = kingpin.Command("covb", "Filter input based on coverage in vcf from same tumor (B).")
	minb     = covb.Flag("min-cov-b", "Minimum coverage in B.").Default("15").Int()
	maxaltb  = covb.Flag("min-alt-b", "Maximum number of alternative alleles in B.").Default("0").Int()
	maxpropb = covb.Flag("min-prop-b", "Maximum proportion of alternative reads in B.").Default("0.0").Float()

	nab      = kingpinCommand("nab", "Filter based on coverage in normal tissue vcf.")
	maxn     = nab.Flag("max-cov-n", "Maximum coverage in normal vcf.").Default("5").Int()
	maxaltn  = nab.Flag("max-alt-n", "Maximum number of alternative alleles in normal vcf.").Default("15").Int()
	maxpropn = nab.Flag("max-prop-n", "Maximum proportion of alternative reads in normal vcf.").Default("0.3").Float()
)

func main() {
	start = time.Now()
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
