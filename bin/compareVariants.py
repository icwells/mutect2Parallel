'''This script will compare variants from different filtering pipelines.'''

import os
import re
import pysam
from shutil import copy
from datetime import datetime
from argparse import ArgumentParser
from sys import stderr
from glob import glob
from commonUtil import *
from unixpath import *

A = re.compile(r"AfiltcovBNAB.*different\.vcf")
B = re.compile(r"BfiltcovBNAB.*different\.vcf")
C = re.compile(r"filtcovBNABU.*common\.vcf")

def comparePipelines(outdir, samples):
	# Calls bcftools isec on each pair of vcfs and writes summary file
	print("\tComparing samples...")
	log = outdir + "comparisonSummary.csv"
	with open(log, "w") as out:
		# Initialize summary file and write header
		out.write("ID,Comparison,Mutect(A),Platypus(B),#PrivateA,#PrivateB,#Common,%Similarity\n")
		for s in samples.keys():
			outpath = outdir + s + "/"
			for t in ["A", "B", "Common"]:
				# Iterate through list to keep fixed order
				a = None
				if "NA" not in samples[s][t] and None not in samples[s][t]:
					a = bcfIsec(outpath + t, samples[s][t])
					if a:
						b = getTotal(outpath + t + "/0001.vcf")
						c = getTotal(outpath + t + "/0002.vcf")
						try:
							sim = c/(a+b+c)
						except ZeroDivisionError:
							sim = 0.0
					if not a:
						a = "NA"
						b = "NA"
						c = "NA"
						sim = 0.0
					out.write(("{},{},{},{},{},{},{},{:.2%}\n").format(s, t,samples[s][t][0], samples[s][t][1], a, b, c, sim))

def comparisonManifest(infile):
	# Reads in dict of vcfs to compare
	samples = {}
	first = True
	with open(infile, "r") as f:
		for line in f:
			if first == False:
				spl = line.strip().split(",")
				if len(spl) == 4:
					s = spl[0]
					if s not in samples.keys():
						samples[s] = {}
					# samples{ID: {type: [vcf1, vcf2]}}
					samples[s][spl[1]] = spl[2:]
			else:
				first = False
	return samples

#-----------------------------------------------------------------------------

def mergeSamples(outfile, mut, plat):
	# Writes samples to manifest file
	print("\tWriting new manifest file...")
	with open(outfile, "w") as out:
		out.write("ID,Sample,Mutect,Platypus\n")
		for i in mut.keys():
			# Write partial lines so they can be editted later
			if i not in plat.keys():
				plat[i] = {}
			for k in ["a", "b", "c"]:
				if k not in plat[i].keys() or plat[i][k] is None:
					plat[i][k] = ""
				if mut[i][k] is None:
					mut[i][k] = ""
			out.write(",".join([i, "A", mut[i]["a"], plat[i]["a"]]) + "\n")
			out.write(",".join([i, "B", mut[i]["b"], plat[i]["b"]]) + "\n")
			out.write(",".join([i, "Common", mut[i]["c"], plat[i]["c"]]) + "\n")

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

def platypusPaths(p):
	# Returns paths from vcfdict
	count = 0
	paths = {}
	with open(p + "vcfdict.csv") as f:
		for line in f:
			spl = line.strip().split(",")
			if count == 3:
				break
			elif A.match(spl[0]):
				paths["a"] = p + spl[1]
				count += 1
			elif B.match(spl[0]):
				paths["b"] = p + spl[1]
				count += 1
			elif C.match(spl[0]):
				paths["c"] = p + spl[1]
				count += 1
	return paths

def getPlatypusOutput(path, outdir = None, contigs = None):
	# Returns dict of platypus output
	plat = {}
	paths = glob(path + "*/")
	print("\tGetting platypus output...")
	for p in paths:
		sdir = None
		p = checkDir(p)
		sample = getParent(p)
		paths = platypusPaths(p)
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

def checkVCF(path, com = False):
	# Determines if file exists and returns filename/NA
	ret = None
	if os.path.isdir(path):
		if com == True:
			path += "common_nab.vcf"
		else:
			path += "/0000.vcf"
		ret = tabix(path, force = True)
	return ret

def getMutectOutput(path):
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
		# Get path names with full sample name
		pa = checkVCF("{}{}".format(p, "A_nab"))
		if pa:
			mut[sample]["a"] = bcfSort(pa)
		pb = checkVCF("{}{}".format(p, "B_nab"))
		if pb:
			mut[sample]["b"] = bcfSort(pb)
		common = checkVCF(p, True)
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
	return args

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will compare variants from different filtering pipelines.")
	parser.add_argument("-c", help = "Copy target platypus data to this directory.")
	parser.add_argument("-v", help = "Path to uncompressed vcf header (Copies contig information to platypus vcf headers).")
	parser.add_argument("-m", help = "Path to mutect2parallel parent output directory.")
	parser.add_argument("-p", help = "Path to platypus-based parent output directory.")
	parser.add_argument("-o", 
help = "Path to output manifest if using -m and -p. Path to output directory if using -i.")
	parser.add_argument("-i", help = "Path to input manifest for comparison.")
	args = checkArgs(parser.parse_args())
	if args.m and args.p:
		print("\n\tGetting new manifest for comparison...")
		contigs = None
		if args.v:
			contigs = getContigs(args.v)
		mutect = getMutectOutput(args.m)
		plat = getPlatypusOutput(args.p, args.c, contigs)
		mergeSamples(args.o, mutect, plat)
	elif args.i:
		# Run comparison
		print("\n\tComparing output from each pipeline...")
		samples = comparisonManifest(args.i)
		comparePipelines(args.o, samples)
	print(("\tFinished. Runtime: {}\n").format(datetime.now()-start))

if __name__ == "__main__":
	main()
