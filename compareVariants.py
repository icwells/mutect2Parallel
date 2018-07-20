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

def getTotal(vcf):
	# Returns total number of content lines from vcf
	count = 0
	if os.path.isfile(vcf):
		with open(vcf, "r") as f:
			for line in f:
				if line[0] != "#":
					count += 1
	return count

def bcfIsec(outpath, vcfs):
	# Calls bcftools to get intersecting rows and summarizes output
	cmd = ("bcftools isec {} {} -p {}").format(vcfs[0], vcfs[1], outpath)
	try:
		bcf = Popen(split(cmd))
		bcf.communicate()
	except:
		print(("\t[Error] Could not call bcftools isec with {}").format(cmd))
		return None
	# Number of unique variants to each sample and number of shared
	a = getTotal(outpath + "/0000.vcf")
	b = getTotal(outpath + "/0001.vcf")
	c = getTotal(outpath + "/0002.vcf")
	# Get percentage of similarity
	try:
		sim = c/(a+b+c)
	except ZeroDivisionError:
		sim = 0.0
	sa = vcfs[0][vcfs[0].rfind("/")+1:vcfs[0].rfind(".")]
	sb = vcfs[1][vcfs[1].rfind("/")+1:vcfs[1].rfind(".")]
	with open(log, "a") as output:
		output.write(("{},{},{},{},{},{}\n").format(sa, sb, a, b, c, sim))
	return 1

def compareVCFS(outpath, vcfs):
	# Calls gatk and pyvcf to filter and compare given pair of
	for i in range(len(vcfs)):
		if vcfs[i] and os.path.isfile(vcfs[i]):
			# Make sure files are bgzipped
			vcfs[i] = bgzip(vcfs[i])
		elif ".gz" not in vcfs[i] and os.path.isfile(vcfs[i] + ".gz"):
			vcfs[i] += ".gz"
		else:
			printError(("Cannot find {}").format(vcfs[i]))
			return None
	res = bcfIsec(outpath,  vcfs)
	return res

def comparePipelines(outdir, samples):
	# Calls bcftools isec on each pair of vcfs and writes summary file
	print("\tComparing samples...")
	log = outdir + "comparisonSummary.csv"
	with open(log, "w") as out:
		# Initialize summary file and write header
		out.write("ID,Comparison,SampleA,SampleB,#PrivateA,#PrivateB,#Common,%Similarity\n")
		for s in samples.keys():
			outpath = outdir + s + "/"
			for t in ["A", "B", "Common"]:
				# Iterate through list to keep fixed order
				res = compareVCFS(outpath + t, samples[s][t])
				if res:
					out.write(("{},{},{}\n").format(s, t, res))


def comparisonManifest(infile):
	# Reads in dict of vcfs to compare
	samples = {}
	first = True
	with open(infile, "r") as f:
		for line in f:
			if first == False:
				spl = line.split(",")
				s = spl[0]
				if s not in samples.keys():
					samples[s] = {}
				# samples{ID: {type: [vcf1, vcf2]}}
				samples[s][spl[1]] = spl[2:]
			else:
				first = False
	return samples

#-----------------------------------------------------------------------------

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

#-----------------------------------------------------------------------------

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
	else:
		args.o = checkDir(args.o, make = True)
	return args

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will compare variants from different filtering pipelines.")
	parser.add_argument("-i", 
help = "Path to input manifest (mutect input for manifest generation or generated manifest for comparison).")
	parser.add_argument("-m", help = "Path to mutect2parallel parent output directory.")
	parser.add_argument("-p", help = "Path to platypus-based parent output directory.")
	parser.add_argument("-o", 
help = "Path to output manifest if using -m and -p. Path to output directory if using -i.")
	args = checkArgs(parser.parse_args())
	if args.m and args.p:
		print("\n\tGetting new manifest for comparison...")
		samples = readManifest(args.i)
		mutect = getMutectOutput(args.m, samples)
		plat = getPlatypusOutput(args.p)
		mergeSamples(args.o, mutect, plat)
	elif args.i:
		# Run comparison
		print("\n\tComparing output from each pipeline...")
		samples = comparisonManifest(args.i)
		comparePipelines(args.o, samples)
	print(("\tFinished. Runtime: {}\n").format(datetime.now()-start))

if __name__ == "__main__":
	main()
