'''This script will compare variants from different filtering pipelines.'''

import os
import re
import pysam
from shutil import copy
from datetime import datetime
from argparse import ArgumentParser
from sys import stderr
from glob import glob
from compareNormals import VCFcomparison, compareSamples
from commonUtil import *
from unixpath import *

A = re.compile(r"AfiltcovBNAB.*different\.vcf")
B = re.compile(r"BfiltcovBNAB.*different\.vcf")
C = re.compile(r"filtcovBNABU.*common\.vcf")

def comparePipelines(samples):
	# Calls bcftools isec on each pair of vcfs and writes summary file
	print("\tComparing samples...")
	for s in samples.keys():
		for t in ["A", "B", "Common"]:
			# Iterate through list to keep fixed order
			if t in samples[s].keys():
				_ = compareSamples(samples[s][t])

def comparisonManifest(infile, outdir):
	# Reads in dict of vcfs to compare
	samples = {}
	first = True
	log = outdir + "comparisonSummary.csv"
	with open(log, "w") as out:
		# Initialize summary file and write header
		out.write("ID,Comparison,Mutect(A),Platypus(B),#PrivateA,#PrivateB,#Common,%Similarity\n")
	# Get input
	with open(infile, "r") as f:
		for line in f:
			if first == False:
				spl = line.strip().split(",")
				if len(spl) == 4:
					s = spl[0]
					if s not in samples.keys():
						samples[s] = {}
					# Store comparison object byt ID and type
					samples[s][spl[1]] = VCFcomparison(spl[1], spl[2], spl[3], log, outdir, sample = spl[0]) 
			else:
				first = False
	return samples

#-----------------------------------------------------------------------------

def isecSamples(out, sample, m, p):
	# Writes a, b, and common pairs to file
	for k in ["a", "b", "c"]:
		if k not in p.keys() or p[k] is None:
			p[k] = ""
		if m[k] is None:
			m[k] = ""
	out.write(",".join([sample, "A", m["a"], p["a"]]) + "\n")
	out.write(",".join([sample, "B", m["b"], p["b"]]) + "\n")
	out.write(",".join([sample, "Common", m["c"], p["c"]]) + "\n")

def rawSamples(out, sample, m, p):
	# Writes sample:common pairs to file
	common = ""
	if "c" in p.keys():
		common = p["c"]
	out.write((",").join([sample, "A", m["a"], common]) + "\n")
	out.write((",").join([sample, "B", m["b"], common]) + "\n")

def mergeSamples(outfile, mut, plat, ext):
	# Writes samples to manifest file
	print("\tWriting new manifest file...")
	with open(outfile, "w") as out:
		out.write("ID,Sample,Mutect,Platypus\n")
		for i in mut.keys():
			# Write partial lines so they can be editted later
			if i not in plat.keys():
				plat[i] = {}
			if ext:
				isecSamples(out, i, mut[i], plat[i])
			else:
				rawSamples(out, i, mut[i], plat[i])

def reheader(contigs, infile, outdir = None):
	# Inserts contig lines into vcf header and returns outfile name
	ins = True
	ids = []
	if outdir:
		outfile = outdir + infile[infile.rfind("/")+1:]
	else:
		outfile = infile.replace(".vcf", "reheadered.vcf")
	# Get list of ids in file and sort
	with open(infile, "r") as f:
		for line in f:
			if line[0] != "#":
				n = line.split("\t")[0]
				if n not in ids:
					ids.append(n)
	with open(outfile, "w") as out:
		with open(infile, "r") as f:
			for line in f:
				if "##FILTER=" in line and ins == True:
					# Insert contigs before first filter line
					for i in ids:
						if i in contigs.keys():
							out.write(contigs[i])
						else:
							print(("\t[Warning] {} not in example vcf header.").format(i), file=stderr)
					ins = False
				out.write(line)
	return outfile

def platypusPaths(p, ext):
	# Returns paths from vcfdict
	count = 0
	paths = {}
	with open(p + "vcfdict.csv") as f:
		for line in f:
			spl = line.strip().split(",")
			if count == 3:
				break
			if A.match(spl[0]) and ext:
				paths["a"] = p + spl[1]
				count += 1
			elif B.match(spl[0]) and ext:
				paths["b"] = p + spl[1]
				count += 1
			elif C.match(spl[0]):
				paths["c"] = p + spl[1]
				if ext:
					count += 1
				else:
					break
	return paths

def getPlatypusOutput(path, outdir, contigs, ext):
	# Returns dict of platypus output
	plat = {}
	paths = glob(path + "*/")
	print("\tGetting platypus output...")
	for p in paths:
		sdir = None
		p = checkDir(p)
		sample = getParent(p)
		paths = platypusPaths(p, ext)
		plat[sample] = {}
		if outdir:
			# Copy to new location
			sdir = checkDir(outdir + sample, True)
		for i in paths.keys():
			if contigs:
				vcf = reheader(contigs, paths[i], sdir)
			elif outdir:
				vcf = copy(paths[i], sdir)
			bcf = bcfSort(tabix(vcf))
			if bcf == None:
				bcf = ""
			plat[sample][i] = bcf 
	return plat

def checkVCF(f):
	# Determines if file exists and returns filename/NA
	ret = tabix(f, force = True)
	return ret

def getNames(path, ext):
	# Returns vcfs names
	c = None
	if ext:
		a = ("{}A_{}/0000.vcf").format(path, ext)
		b = ("{}B_{}/0000.vcf").format(path, ext)
		c = ("{}common_{}.vcf").format(path, ext)
	else:
		a = ("{}A.vcf").format(path)
		b = ("{}B.vcf").format(path)
	return a, b, c

def getMutectOutput(path, ext):
	# Returns dict of fitered mutect output
	mut = {}
	paths = glob(path + "*/")
	print("\tGetting mutect output...")
	for p in paths:
		p = checkDir(p)
		full = getParent(p)
		# Drop leading number from name
		sample = full[full.find("_")+1:]
		mut[sample] = {}
		mut[sample]["a"] = "NA"
		mut[sample]["b"] = "NA"
		mut[sample]["c"] = "NA"
		# Get path names and check each file
		a, b, c = getNames(p, ext)
		pa = checkVCF(a)
		if pa:
			mut[sample]["a"] = bcfSort(pa)
		pb = checkVCF(b)
		if pb:
			mut[sample]["b"] = bcfSort(pb)
		if c:
			common = checkVCF(c)
			if common:
				mut[sample]["c"] = bcfSort(common)
	return mut

#-----------------------------------------------------------------------------

def getContigs(infile):
	# Reads contig info from vcf header
	contigs = {}
	print("\tGetting contigs from vcf header...")
	with open(infile, "r") as f:
		for line in f:
			if line[0] != "#":
				break
			if "##contig=" in line:
				spl = line.split(",")
				name = spl[0][spl[0].rfind("=")+1:]
				contigs[name] = line
	return contigs

def getExt(raw, unfiltered):
	# Returns extension for type of comparison
	if raw == True:
		return None
	elif unfiltered == True:
		return "unfiltered"
	return "nab"

def checkArgs(args):
	# Checks for errors in arguments
	if args.m or args.p:
		if not args.m or not args.p:
			# Raise error if either is missing
			print("\n\t[Error] Please specify both -m and -p if generating a manifest. Exiting.\n", file=stderr)
			quit()
		args.m = checkDir(args.m)
		args.p = checkDir(args.p)
		if not args.o or os.path.isdir(args.o):
			print("\n\t[Error] Please specify an output file. Exiting.\n", file=stderr)
			quit()
		if args.c:
			args.c = checkDir(args.c, True)
		if args.v:
			checkFile(args.v)
	else:
		checkFile(args.i)
		args.o = checkDir(args.o, True)
	return args, getExt(args.raw, args.unfiltered)

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will compare variants from different filtering pipelines.")
	parser.add_argument("--raw", action = "store_true", default = False,
help = "Compare Raw mutect2 output to common variants from platypus output (compares fully filtered output by default).")
	parser.add_argument("--unfiltered", action = "store_true", default = False,
help = "Compare 'unfiltered' somatic variants (compares fully filtered output by default).")
	parser.add_argument("-c", help = "Copy target platypus data to this directory.")
	parser.add_argument("-v", help = "Path to uncompressed vcf header (Copies contig information to platypus vcf headers).")
	parser.add_argument("-m", help = "Path to mutect2parallel parent output directory.")
	parser.add_argument("-p", help = "Path to platypus-based parent output directory.")
	parser.add_argument("-o", 
help = "Path to output manifest if using -m and -p. Path to output directory if using -i.")
	parser.add_argument("-i", help = "Path to input manifest for comparison.")
	args, ext = checkArgs(parser.parse_args())
	if args.m and args.p:
		print("\n\tGetting new manifest for comparison...")
		contigs = None
		if args.v:
			contigs = getContigs(args.v)
		mutect = getMutectOutput(args.m, ext)
		plat = getPlatypusOutput(args.p, args.c, contigs, ext)
		mergeSamples(args.o, mutect, plat, ext)
	elif args.i:
		# Run comparison
		print("\n\tComparing output from each pipeline...")
		samples = comparisonManifest(args.i, args.o)
		comparePipelines(samples)
	print(("\tFinished. Runtime: {}\n").format(datetime.now()-start))

if __name__ == "__main__":
	main()
