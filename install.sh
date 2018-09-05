#!/bin/bash

##############################################################################
# This script will install scripts for the compOncDB package.
# 
# Required programs:	Go 1.7+
##############################################################################

IO="github.com/icwells/go-tools/iotools"
KP="gopkg.in/alecthomas/kingpin.v2"
MAIN="heterAnalyzer"
PR="github.com/Songmu/prompter"
SA="github.com/icwells/go-tools/strarray"

# Get install location
SYS=$(ls $GOPATH/pkg | head -1)
PDIR=$GOPATH/pkg/$SYS

echo ""
echo "Preparing compOncDB package..."
echo "GOPATH identified as $GOPATH"
echo ""

# Get dependencies
for I in $IO $KP $PR $SA; do
	if [ ! -e "$PDIR/$I.a" ]; then
		echo "Installing $I..."
		go get -u $I
		echo ""
	fi
done

# lineageSimulator 
echo "Building main..."
go build -o bin/$MAIN src/$MAIN/*.go

echo "Finished"
echo ""
