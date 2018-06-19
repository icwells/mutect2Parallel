'''This script will copy all csv summaries from filterVCFs into one summary file'''

import os
from argparse import ArgumentParser
from datetime import datetime
from glob import glob

def getCSV(csv):
	# Returns list of summary files
	with open(csv, "r") as f:
		s = f.readlines()
	# Return file data without header
	return s[1:]

def summarize(indir, outfile):
	# Copies contents of each summary csv in each subdirectory
	paths = glob(conf["outpath"] + "*")
	with open(outfile, "w") as out:
		out.write("Sample,Type,SampleA,SampleB,#PrivateA,#PrivateB,#Common,Similarity\n")
		for p in paths:
			# Iterate through each subdirectory
			if p[:-1] != "/":
				p += "/"
			# Get sample name from path
			sample = p.split("/")[-1]
			csv = p + sample + ".csv"
			if not os.path.isfile(csv):
				print(("\t[Error] Cannot find {}. Skipping.").format(csv))
			else:
				s = getCSV(csv)
				if s:
					for i in s:
						# Write sample name and original line
						out.write(("{},{}\n").format(Sample, i))

def main():
	start = datetime.now()
	parser = ArgumentParser(
"This script will copy all csv summaries from filterVCFs into one summary file")
	parser.add_argument("-i", help = "Path to mutect2Parallel output directory.")
	parser.add_argument("-o", help = "Path to output csv.")
	args = parser.parse_args()
	if not args.i or not args.o:
		print("\n\t[Error] Please specify input directory and output file. Exiting.\n")
		quit()
	if args.i[-1] != "/":
		args.i += "/"
	summarize(args.i, args.o)
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-start))	

if __name__ == "__main__":
	main()
