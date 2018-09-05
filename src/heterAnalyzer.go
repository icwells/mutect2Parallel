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

var(
	app    = kingpin.New("heterAnalyzer", "Applies custom filters to vcf files.")
	infile = kingpin.Flag("vcf", "Path to input vcf.").Required().Short('v').Default("nil").String()
	interval = kingpin.Flag("interval", "Path to covB.sh/nab.sh output.").Required().Short('i').Default("nil").String()
	outfile = kingpin.Flag("outfile", "Path to output vcf.").Required().Short('o').Default("nil").String()

	covb = kingpin.Command("covb", "Filter input based on coverage in vcf from same tumor (B).")
	minb = covb.Flag("min-cov-b", "Minimum coverage in B.").Default("15").Int()
	minaltb = covb.Flag("min-alt-b", "Minimum number of alternative alleles in B.").Default("0").Int()
	minpropb = covb.Flag("min-prop-b", "Minimum proportion of alternative reads in B.").Default("0.0").Float()

	nab = kingpinCommand("nab", "Filter based on coverage in normal tissue vcf.")
	maxcovn = nab.Flag("max-cov-n", "Maximum coverage in normal vcf.").Default("5").Int()
	maxaltn = nab.Flag("max-alt-n", "Maximum number of alternative alleles in normal vcf.").Default("15").Int()
	maxpropn = nab.Flag("max-prop-n", "Maximum proportion of alternative reads in normal vcf.").Default("0.3").Float()
)

func main() {
	start = time.Now()
	var vcf VCF
	// Get command to append to vcf
	cmd := strings.Join(os.Args(), " ")
	switch kingpin.Parse() {
		case covb.FullCommand():
			vcf.readVCF(*infile, cmd)
			filter
		case nab.FullCommand():
			vcf := readVCF(*infile)
	}
	fmt.Printf("\n\tFinished. Runtime: %s\n\n", time.Since(start))
}
