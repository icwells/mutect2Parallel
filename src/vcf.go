// This script defines a struct for storing a vcf file

package main

import (
	"fmt"
	"github.com/icwells/go-tools/iotools"
	"os"
	"strings"
)

type VCF struct {
	Header []string
	Regions map[string]map[int][]string
}

func (v *VCF) readVCF(infile, cmd string) {
	// Reads vcf file into struct
	v.Regions = make(map[string]map[int][]string)
	fmt.Printf("\n\tReading vcf: %s\n", infile)
	f := iotools.OpenFile(infile)
	defer f.Close()
	for input.Scan() {
		line := string(input.Text())
		if strings.Contains(line, "##") == true {
			// Append raw header lines
			v.Header = append(v.Header, line)
		} else if strings.Contains(line, "#") == true {
			// Append command prior to column header line
			v.Header = append(v.Header, fmt.Sprintf("##%s\n", cmd))
			v.Header = append(v.Header, line)
		} else {
			s := strings.Split(line, "\t")
			
		}
	}
}
