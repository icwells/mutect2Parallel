#!/bin/bash

##############################################################################
# This script will install scripts for mutect2Parallel.
# 
# Required programs:	Go 1.7+
##############################################################################

FI="filterNAB"
IO="github.com/icwells/go-tools/iotools"
KP="gopkg.in/alecthomas/kingpin.v2"
SA="github.com/icwells/go-tools/strarray"

# Get install location
SYS=$(ls $GOPATH/pkg | head -1)
PDIR=$GOPATH/pkg/$SYS

echo ""
echo "Preparing filterNAB package..."
echo "GOPATH identified as $GOPATH"
echo ""

# Get dependencies
for I in $IO $KP $SA; do
	if [ ! -e "$PDIR/$I.a" ]; then
		echo "Installing $I..."
		go get -u $I
		echo ""
	fi
done

# lineageSimulator 
echo "Building main..."
go build -o bin/$FI src/*.go

echo "Finished"
echo ""
