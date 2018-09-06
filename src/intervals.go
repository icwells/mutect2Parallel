// This script defines a struct for storing a covb/nan output

package main

import (
	"fmt"
	"github.com/icwells/go-tools/iotools"
	"strconv"
	"strings"
)

type Params struct {
	nab     bool
	minb    int
	maxaltb int
	maxprob float64
	minn    int
	maxaltn int
	maxpron float64
}

func (p *Params) setParameters(nab bool, minb, maxaltb, minn, maxaltn int, maxprob, maxpron float64) {
	// Stores parameters for filtering
	p.nab = nab
	p.minb = minb
	p.maxaltb = maxaltb
	p.maxprob = maxprob
	p.minn = minn
	p.maxaltn = maxaltn
	p.maxpron = maxpron
}

type Region struct {
	Ref      string
	Alt      string
	RefReads int
	AltReads int
}

func (r *Region) setRegion(ref, alt, rreads, areads string) {
	// Converts and stores region values
	r.Ref = ref
	r.Alt = alt
	rr, _ := strconv.ParseInt(rreads, 10, 64)
	ar, _ := strconv.ParseInt(areads, 10, 64)
	r.RefReads = int(rr)
	r.AltReads = int(ar)
}

type Intervals struct {
	Settings Params
	Regions  map[string]map[int]Region
}

func (i *Intervals) setIntervals(nab bool, minb, maxaltb, minn, maxaltn int, maxprob, maxpron float64) {
	// Initializes interval values
	i.Regions = make(map[string]map[int]Region)
	i.Settings.setParameters(nab, minb, maxaltb, minn, maxaltn, maxprob, maxpron)
}

func (i *Intervals) loadIntervals(infile string) {
	// Reads interval file into struct
	fmt.Printf("\n\tReading interval file: %s\n", infile)
	f := iotools.OpenFile(infile)
	defer f.Close()
	input := iotools.GetScanner(f)
	for input.Scan() {
		line := string(input.Text())
		s := strings.Split(line, "\t")
		_, exists := i.Regions[s[0]]
		if exists == false {
			// Initialize new map for each chromosome
			i.Regions[s[0]] = make(map[int]Region)
		}
		// Use coordinates as unique map key to store region values
		var reg Region
		reg.setRegion(s[2], s[3], s[4], s[5])
		loc, _ := strconv.ParseInt(s[1], 10, 64)
		i.Regions[s[0]][int(loc)] = reg
	}
}
