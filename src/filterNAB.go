// This script will filter target cancer vcf files in comparison to a normal vcf and an unfiltered paired vcf

package main

import (
	"fmt"
	"github.com/icwells/go-tools/iotools"
	"github.com/icwells/go-tools/strarray"
	"gopkg.in/alecthomas/kingpin.v2"
	//"os"
	//"os/exec"
	"strings"
	"time"
)

func filterA(outfile, infile string, variants map[string][]string) {
	// Writes variants from infile to outfile if they re present in variants
	out := iotools.CreateFile(outfile)
	defer out.CLose()
	f := iotools.OpenFile(infile)
	defer f.Close()
	input := iotools.GetScanner(f)
	for input.Scan() {
		if line[0] != '#' {
			err := out.WriteString(line)
			iotools.CheckError("Writing line to vcf", err, 0)
		} else {
			s := strings.Split(line, "\t")
			l, pass := variants[s[0]]
			if pass == true || strarray.InSliceStr(l, s[1]) == true {
				err := out.WriteString(line)
				iotools.CheckError("Writing line to vcf", err, 0)
			}
		}
	}
}

func sortVariants(conf Config, a, b, c VCF) map[string][]string {
	// Applys setting in variants to approriate vcfs and removes variants in a which do not pass
	ret := make(map[string][]string)

	return ret
}

func main() {
	start := time.Now()
	var (
		_       = kingpin.New("filterNAB", "This script will filter target cancer vcf files in comparison to a normal vcf and an unfiltered paired vcf.")
		c       = kingpin.Flag("config", "PAth to config file.").Short('c').Required.String()
		vcfa    = kingpin.Flag("vcfA", "VCF file to filter.").Short('a').Required().String()
		vcfb    = kingpin.Flag("vcfB", "Unfiltered VCF from the same tumor.").Short('b').Required().String()
		vcfn    = kingpin.Flag("normal", "Normal vcf from the same tumor.").Short('n').Required().String()
		outfile = kingpin.Flag("outfile", "Output vcf file.").Short('o').Required().String()
	)
	kingpin.Parse()
	var conf Config
	var a, b, n VCF
	// Get config and vcfs
	conf.setConfig(*c)
	a.LoadVCF(*vcfa)
	b.LoadVCF(*vcfb)
	n.LoadVCF(*vcfn)
	// Filter and write output
	variants = sortVariants(conf, a, b, n)
	filterA(*outfile, *vcfa, variants)
	fmt.Printf("\n\tFinished. Runtime: %s\n\n", time.Since(start))
}
