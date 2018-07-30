'''This script will compare variants from different filtering pipelines.'''

import os
import re
import pysam
from shutil import copy
from datetime import datetime
from argparse import ArgumentParser
from glob import glob
from shlex import split
from subprocess import Popen
from commonUtil import *
from bamUtils import *
from unixpath import *

A = re.compile(r"AfiltcovBNAB.*different\.vcf")
B = re.compile(r"BfiltcovBNAB.*different\.vcf")
C = re.compile(r"filtcovBNABU.*common\.vcf")

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
	for i in range(len(vcfs)):
		# MAke sure there is an up-to-date index file
		vcfs[i] = tabix(vcfs[i], True)
	if None in vcfs:
		return None
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
	return ("{},{},{},{},{},{}\n").format(sa, sb, a, b, c, sim)

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
				res = bcfIsec(outpath + t, samples[s][t])
				if res:
					out.write(("{},{},{}").format(s, t, res))

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
				if k not in plat[i].keys():
					plat[i][k] = ""
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
							print(("\t[Warning] {} not in example vcf header.").format(i))
					ins = False
				out.write(line)
	return outfile

def getPlatypusOutput(path, outdir = None, contigs = None):
	# Returns dict of platypus output
	plat = {}
	paths = glob(path + "*/")
	print("\tGetting platypus output...")
	for p in paths:
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
				else:
					vcf = copy(paths[i], sdir)
				bcf = bcfSort(tabix(vcf))
				if bcf == None:
					bcf = ""
				plat[sample][i] = bcf 
		else:
			for i in paths.keys():
				if contigs:
					vcf = reheader(contigs, paths[i], sdir)
				else:
					vcf = paths[i]	
				bcf = bcfSort(tabix(vcf))
				if bcf == None:
					bcf = ""
				plat[sample][i] = bcf 
	return plat

def getMutectOutput(path, samples):
	# Returns dict of fitered mutect output
	mut = {}
	paths = glob(path + "*/")
	print("\tGetting mutect output...")
	for p in paths:
		p = checkDir(p)
		full = getParent(p)
		# Drop leading number from name
		sample = full[full.find("_")+1:]
		if sample in samples.keys():
			com = []
			a = samples[sample]["a"]
			b = samples[sample]["b"]
			if len(a) > 1 and len(b) > 1:
				# Get path names with full sample name
				pa = checkDir("{}{}_{}".format(p, full, a))
				pb = checkDir("{}{}_{}".format(p, full, b))
				# Sort and merge common vcfs
				for i in [pa, pb]:
					srt = bcfSort(tabix(i + "0003.vcf"))
					com.append(srt)
				if None not in com:
					common = bcfMerge(p, com)
				if common != None and os.path.isfile(common):
					# Save private filtered A and B, and common
					mut[sample] = {}
					mut[sample]["c"] = common
					mut[sample]["a"] = bcfSort(tabix(pa + "0000.vcf"))
					mut[sample]["b"] = bcfSort(tabix(pb + "0000.vcf"))
				else:
					print(("\t[Error] calling bcftools on {}.").format(sample))
		else:
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
	checkFile(args.i)
	if args.m or args.p:
		if not args.m or not args.p:
			# Raise error if either is missing
			print("\n\t[Error] Please specify both -m and -p if generating a manifest. Exiting.\n")
			quit()
		args.m = checkDir(args.m)
		args.p = checkDir(args.p)
		if not args.o or os.path.isdir(args.o):
			print("\n\t[Error] Please specify an output file. Exiting.\n")
			quit()
		if args.c:
			args.c = checkDir(args.c, True)
		if args.v:
			checkFile(args.v)
	else:
		args.o = checkDir(args.o, True)
	return args

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will compare variants from different filtering pipelines.")
	parser.add_argument("-c", help = "Copy target platypus data to this directory.")
	parser.add_argument("-v", help = "Path to uncompressed vcf header (Copies contig information to platypus vcf headers).")
	parser.add_argument("-i", 
help = "Path to input manifest (mutect input for manifest generation or generated manifest for comparison).")
	parser.add_argument("-m", help = "Path to mutect2parallel parent output directory.")
	parser.add_argument("-p", help = "Path to platypus-based parent output directory.")
	parser.add_argument("-o", 
help = "Path to output manifest if using -m and -p. Path to output directory if using -i.")
	args = checkArgs(parser.parse_args())
	if args.m and args.p:
		print("\n\tGetting new manifest for comparison...")
		contigs = None
		if args.v:
			contigs = getContigs(args.v)
		samples = readManifest(args.i)
		mutect = getMutectOutput(args.m, samples)
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
