'''This script will subset chromosome(s) from a bed file for use as as the active region for Mustect2'''

import os
from datetime import datetime
from argparse import ArgumentParser

def subsetRegions(regions, infile, outfile):
	# Copies entries from given chromosomes to outfile
	print("\n\tSubsetting interval file...")
	with open(outfile, "w") as output:
		with open(infile, "r") as f:
			for line in f:
				splt = line.split("\t")
				for i in regions:
					if i == splt[0]:
						output.write(line)
						break

def getRegions(l):
	# Returns list of chromosomes from input
	if "," in l:
		return l.split(",")
	else:
		return [l]	

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will subset chromosome(s) from a \
bed file for use as as the active region for Mustect2.")
	parser.add_argument("-c",
help = "Chromosome(s) to subset (seperate with commas if there is more than one).")
	parser.add_argument("-i", help = "Path to input file.")
	parser.add_argument("-o", help = "Path to output file.")
	args = parser.parse_args()
	if not os.path.isfile(args.i):
		print(("\n\t[Error] Cannot find {}. Exiting.\n").format(args.i))
		quit()
	regions = getRegions(args.c)
	subsetRegions(regions, args.i, args.o)
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-start))

if __name__ == "__main__":
	main()
