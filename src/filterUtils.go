// This file contains helper functions for filterNAB

package main

import (
	"buifio"
	"fmt"
	"github.com/icwells/go-tools/iotools"
	"os"
	"strings"
)

type Variant struct {
	pos    int
	ref    string
	alt    string
	filter string
}

func (v *Variant) setVariant(s []string) {
	// Reads split line into struct
	v.pos = strconv.ParseInt(s[1], 10, 64)
	v.ref = s[3]
	v.alt = s[4]
	v.filter = s[6]
}

type VCF struct {
	name string
	//header		[]string
	variants map[string][]*Variant
}

func (v *VCF) LoadVCF(infile string) {
	// Reads vcf file into struct
	v.name = iotools.GetFileName(infile)
	f := iotools.OpenFile(infile)
	defer f.Close()
	input := iotools.GetScanner(f)
	for input.Scan() {
		if line[0] != '#' {
			// Assign variants to map of chromosome regions
			s := strings.Split(line, "\t")
			var variant Variant
			variant.setVariant(s)
			v.variants[s[0]] = append(v.variants[s[0]], &variant)
			//} else {
			//	v.header = append(v.header, line)
		}
	}
}

type Config struct {
	qual             int
	min_covA         int
	min_reads_strand int
	min_reads_alt    int
	min_covB         int
	max_altB         int
	max_prop_altB    float64
	max_covN         int
	min_freq_altN    float64
	min_reads_altN   int
}

func (c *Config) setConfig(infile string) {
	// Reads config file and assigns to stuct
	var err error
	if iotools.Exists(infile) == false {
		fmt.Printf("\n\t[Error] Cannot find %s\n", infile)
		os.Exit(2)
	}
	fmt.Println("\n\tReading config file...\n")
	f := iotools.OpenFile(infile)
	defer f.Close()
	input := bufio.NewScanner(f)
	for input.Scan() {
		line := string(input.Text())
		if strings.Contains(line, "#SBATCH") == true || strings.Contains(line, "#PBS") == true {
			// Skip sample batch script
			break
		}
		if strings.Contains(line, "#") == false {
			// Skip commented lines
			s := strings.Split(line, "=")
			k := strings.Trim(s[0], "\n\t ")
			v := strings.Trim(s[1], "\n\t ")
			switch k {
			// Assign values to struct fields
			case "qual":
				c.qual, err = strconv.ParseInt(v, 10, 64)
				_ := iotools.CheckError("qual is not an integer", err, 10)
			case "min_covA":
				c.min_covA, err = strconv.ParseInt(v, 10, 64)
				_ := iotools.CheckError("min_covA is not an integer", err, 10)
			case "min_reads_strand":
				c.min_reads_strand, err = strconv.ParseInt(v, 10, 64)
				_ := iotools.CheckError("min_reads_strand is not an integer", err, 10)
			case "min_reads_alt":
				c.min_reads_alt, err = strconv.ParseInt(v, 10, 64)
				_ := iotools.CheckError("min_reads_alt is not an integer", err, 10)
			case "min_covB":
				c.min_covB, err = strconv.ParseInt(v, 10, 64)
				_ := iotools.CheckError("min_covB is not an integer", err, 10)
			case "max_altB":
				c.max_altB, err = strconv.ParseInt(v, 10, 64)
				_ := iotools.CheckError("max_altB is not an integer", err, 10)
			case "max_prop_altB":
				c.max_prop_altB, err = strconv.ParseFloat(v, 64)
				_ := iotools.CheckError("max_prop_altB is not a float", err, 10)
			case "max_covN":
				c.max_covN, err = strconv.ParseInt(v, 10, 64)
				_ := iotools.CheckError("max_covN is not an integer", err, 10)
			case "min_freq_altN":
				c.min_freq_altN, err = strconv.ParseFloat(v, 64)
				_ := iotools.CheckError("min_freq_altN is not a float", err, 10)
			case "min_reads_altN":
				c.min_reads_altN, err = strconv.ParseInt(v, 10, 64)
				_ := iotools.CheckError("min_reads_altN is not an integer", err, 10)
			}
		}
	}
}
