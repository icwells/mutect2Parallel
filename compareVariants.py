'''This script will compare variants from different filtering pipelines.'''

import os
import re
from datetime import datetime
from argparse import ArgumentParser
from glob import glob
from shlex import split
from subprocess import Popen
from unixpath import *

A = re.compile(r"AfiltcovBNAB.*different.vcf")
B = re.compile(r"BfiltcovBNAB.*different.vcf")
C = re.compile(r"filtcovBNABU.*common.vcf")

def bcfConcat(path, com):
	# Calls bftools concat on input files
	outfile = path + "common.vcf"
	cmd = ("bcftools concat -D -O v -o {} {} {}").format(outfile, com[0], com[1])
	with open(os.devnull, "w") as dn:
		try:
			bc = Popen(split(cmd), stdout = dn, stderr = dn)
			bc.communicate()
		except:
			print(("\t[Error] calling bcftools concat on samples in {}").format(path))
			outfile = None
	return outfile	

def bcfSort(infile):
	# Call bcftools sort
	outfile = infile.replace(".vcf", ".sorted.vcf")
	cmd = ("bcftools sort -O v -o {} {}").format(outfile, infile)
	with open(os.devnull, "w") as dn:
		try:
			bs = Popen(split(cmd), stdout = dn, stderr = dn)
			bs.communicate()
		except:
			print(("\t[Error] calling bcftools sort on {}").format(infile))
			outfile = None
	return outfile

#-----------------------------------------------------------------------------

def mergeSamples(outfile, mut, plat):
	# Writes samples to manifest file
	print("\tWriting new manifest file...")
	with open(outfile, "w") as out:
		out.write("ID,Sample,Mutect,Platypus\n")
			for i in mut.keys():
				if i in plat.keys():
					out.write(",".join([i, "A", mut[i]["a"], plat[i]["a"]]) + "\n")
					out.write(",".join([i, "B", mut[i]["b"], plat[i]["b"]]) + "\n")
					out.write(",".join([i, "Common", mut[i]["c"], plat[i]["c"]]) + "\n")

def platypusPaths(p):
	# Returns paths from vcfdict
	count = 0
	paths = {}
	with open(p + "vcfdict.csv") as f:
		for line in f:
			spl = line.strip().split(",")
			if count == 3:
				break
			elif A.match(spl[0]) == True:
				paths["a"] = spl[1]
				count += 1
			elif B.match(spl[0]) == True:
				paths["b"] = spl[1]
				count += 1
			elif C.match(spl[0]) == True:
				paths["c"] = spl[1]
				count += 1
	return paths	

def getPlatypusOutput(path):
	# Returns dict of platypus output
	plat = {}
	paths = glob(path + "*/")
	print("\tGetting platypus output...")
	for p in paths:
		count = 
		p = checkDir(p)
		sample = getParent(p)
		plat[sample] = platypusPaths(p)
	return plat

def getMutectOutput(path, samples):
	# Returns dict of fitered mutect output
	mut = {}
	paths = glob(path + "*/")
	print("\tGetting mutect output...")
	for p in paths:
		p = checkDir(p)
		sample = getParent(p)
		# Drop leading number from name
		sample = sample[sample.find("_")+1:]
		if sample in samples.keys():
			com = []
			a = samples[sample]["a"]
			b = samples[sample]["b"]
			pa = checkDir(p + a)
			pb = checkDir(p + b)
			# Sort and concatenate common vcfs
			for i in [pa, pb]:
				srt = bcfSort(i + "0003.vcf")
				com.append(srt)
			if None not in com:
				common = bcfConcat(p, com)
			if common != None:
				# Save private filtered A and B, and common
				mut[sample] = {}
				mut[sample]["c"] = common
				mut[sample]["a"] = pa + "0000.vcf"
				mut[sample]["b"] = pb + "0000.vcf"
		print(("\t[Error] {} not in manifest.").format(sample))
	return mut

def readManifest(infile):
	# Returns dict of sample names
	samples = {}
	first = True
	print("\tReading mutect manifest...")
	with open(infile, "r") as f:
		for line in f:
			if line[0] != "#" and line.split():
				if first == True:
					delim = getDelim(line)
					first = False
				spl = line.strip().split(delim)
				s = spl[0][spl[0].find("_")+1:]
				a = getFileName(spl[2])
				b = getFileName(spl[3])
				samples[s] = {"a": a, "b": b}
	return samples

def checkArgs(args):
	# Checks for errors in arguments
	args.i = checkFile(args.i)
	if args.m or args.p:
		if not args.m or not args.p:
			# Raise error if either is missing
			print("\n\t[Error] Please specify both -m and -p if generating a manifest. Exiting.\n")
			quit()
		args.m = checkDir(args.m)
		args.p = checkDir(args.p)
		if not args.o or os.path.isdir(args.o):
			print("\n\t[Error] Please specify and output file. Exiting.\n")
			quit()
		args.o = checkFile(args.o)
	elif:
		args.o = checkDir(args.o, make = True)
	return args

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will compare variants from different filtering pipelines.")
	parser.add_argument("-i", help = "Path to input manifest.")
	parser.add_arguments("-m", help = "Path to mutect2parallel parent output directory.")
	parser.add_argument("-p", help = "Path to platypus-based parent output directory.")
	parser.add_argument("-o", help = "Path to output manifest if using -m and -p. Path to output directory if using -i.")
	args = checkArgs(parser.parser_args())
	if args.m and args.p:
		print("\n\tGetting new manifest for comparison...")
		samples = readManifest(args.i)
		mutect = getMutectOutput(args.m, samples)
		plat = getPlatypusOutput(args.p)
		mergeSamples(args.o, mutect, plat)
	elif args.i:
		# Run comparison
		print("\n\tComparing output from each pipeline...")

if __name__ == "__main__":
	main()
