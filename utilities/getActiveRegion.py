'''This script will subset chromosome(s) from a bed file for use as as the active region for Mutect2'''

import os
from datetime import datetime
from argparse import ArgumentParser

def subsetRegions(regions, s, infile, outfile, outfile2):
	# Copies entries from given chromosomes to outfile
	print("\n\tSubsetting interval file...")
	nomatch = []
	other = set()
	with open(outfile, "w") as output:
		with open(infile, "r") as f:
			for line in f:
				match = False
				splt = line.split("\t")
				for i in regions:
					if i == splt[0]:
						output.write(line)
						match = True
						break
				if match == False:
					if s == True:
						other.add(splt[0])
						nomatch.append(line)
	if s == True:
		other = list(other)
		other.sort()
		print(("\tWriting the following chromosome to {}...\n").format(outfile2))
		print(other)
		with open(outfile2, "w") as out:
			for line in nomatch:
				out.write(line)

def getRegions(l):
	# Returns list of chromosomes from input
	if "," in l:
		return l.split(",")
	else:
		return [l]

def checkArgs(args):
	# Checks arguments for input errors
	if not args.i or not args.o:
		print("\n\t[Error] Please specify input and output files. Exiting.\n", file = sys.stderr)
		quit()
	if not os.path.isfile(args.i):
		print(("\n\t[Error] Cannot find {}. Exiting.\n").format(args.i))
		quit()
	if args.split == True and not args.n:
		print("\n\t[Error] Please specify second output file with -n. Exiting.\n", file = sys.stderr)
		quit()

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will subset chromosome(s) from a \
bed file for use as as the active region for Mustect2.")
	parser.add_argument("--split", action = "store_true", default = False,
help = "Additionally write sequences which do match to seperate output file.")
	parser.add_argument("-c",
help = "Chromosome(s) to subset (seperate with commas if there is more than one).")
	parser.add_argument("-i", help = "Path to input file.")
	parser.add_argument("-o", help = "Path to output file.")
	parser.add_argument("-n", help = "Path to file of non-maached sequences if using --split.")
	args = parser.parse_args()
	checkArgs(args)
	regions = getRegions(args.c)
	subsetRegions(regions, args.split, args.i, args.o, args.n)
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-start))

if __name__ == "__main__":
	main()
