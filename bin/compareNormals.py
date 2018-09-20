'''This script will call bcftools isec to identify any mis-labled samples'''

import os
from datetime import datetime
from argparse import ArgumentParser
from shutil import copy
from itertools import combinations
from multiprocessing import Pool, cpu_count
from functools import partial
from unixpath import *
from commonUtil import *

def identifySample(norm, outdir, log, vcfs):
	# Calls bcftools isec for each set of normals vs vcf
	v = getFileName(vcfs[0])
	n = getFileName(vcfs[1])
	if norm == True:
		outpath = ("{}isecNormals/{}_{}").format(outdir, v, n)
	else:
		outpath = ("{}isecSample/{}_{}").format(outdir, v, n)
	a = bcfIsec(outpath, vcfs)
	if a is not None:
		b = getTotal(outpath + "/0001.vcf")
		c = getTotal(outpath + "/0002.vcf")
		try:
			sim = c/(a+b+c)
		except ZeroDivisionError:
			sim = 0.0
		with open(log, "a") as out:
			out.write(("{},{},{},{},{},{:.2%}\n").format(v, n, a, b, c, sim))
		return [True, v, n]
	else:
		return [False, v, n]

def allSamplePairs(outdir, normals, a, b):
	# Returns all pairs for a:normal and b:normal
	vcfs = []
	log = outdir + "allSamplesComparison.csv"
	with open(log, "w") as out:
		# Initialize log file
		out.write("SampleA, SampleB, PrivateA, PrivateB, Common, %Similarity\n")
	for t in [a, b]:
		for i in t:
			for j in normals:
				# Append each a/b to normal pair
				vcfs.append([i, j])
	return vcfs, log

def getSamplePairs(outdir, normals, vcf = None):
	# Returns pairs of samples to compare
	vcfs = []
	if vcf:
		log = outdir + getFileName(vcf) + "Comparison.csv"
		for i in normals:
			# Pair input vcf with each normal vcf
			vcfs.append([vcf, i])
	else:
		log = outdir + "normalsComparison.csv"
		c = combinations(normals, 2)
		# Convert to list
		for i in c:
			v = [i[0], i[1]]
			vcfs.append(v)
	with open(log, "w") as out:
		# Initialize log file
		out.write("SampleA, SampleB, PrivateA, PrivateB, Common, %Similarity\n")
	return vcfs, log

def checkVCF(line, outpath, stem):
	# Returns name of existing output file and copies if necessary
	outfile = None
	if stem not in line:
		line = line.replace(".N.vcf", stem)
	if os.path.isfile(line):
		outfile = outpath + getParent(line) + stem
		if os.path.isfile(outfile + ".gz"):
			# Check for existance of gzipped file
			outfile += ".gz"
		elif not os.path.isfile(outfile):
			# Copy file if it has not already been copied
			copy(line, outfile)
	else:
		print(("\t[Warning] {} not found. Skipping.").format(line))
	return outfile

def getNormals(infile, outdir, allsamples):
	# Reads input manifest and copies files to outdir
	normals = []
	tumora = []
	tumorb = []
	outpath = checkDir(outdir + "vcfs/", True)
	print("\n\tReading normals manifest...")
	with open(infile, "r") as f:
		for line in f:
			line = line.strip()
			if allsamples == False:
				outfile = checkVCF(line, outpath, ".N.vcf")
				if outfile:
					normals.append(outfile)
			else:
				for i in [".N.vcf", ".A.vcf", ".B.vcf"]:
					# Get outfile name with sample as part of file name
					outfile = checkVCF(line, outpath, i)
					if outfile:
						if "N" in i:
							normals.append(outfile)
						elif "A" in i:
							tumora.append(outfile)
						elif "B" in i:
							tumorb.append(outfile)
	return normals, tumora, tumorb

def fatalError(msg):
	# Prints meassage and exits
	print(("\n\t[Error] {}. Exiting.\n").format(msg))
	quit()

def checkArgs(args):
	# Makes sure required input are given and input files exist
	norm = False
	if not args.o:
		fatalError("Output directory required")
	elif not args.m:
		fatalError("Manifest of normals required")
	if args.i:
		if os.path.isfile(args.i + ".gz"):
			args.i = args.i + ".gz"
		else:
			checkFile(args.i)
	else:
		norm = True
	checkFile(args.m)
	args.o = checkDir(args.o, True)
	if args.t > cpu_count():
		args.t = cpu_count()
	return args, norm

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will call bcftools isec to compare input samples. \
Make sure platypus is loaded in a module or in your PATH if supplying an input bam file.")
	parser.add_argument("--allsamples", action = "store_true", default = False,
help = "Compares each tumor vcf to all normal vcfs.")
	parser.add_argument("-t", type = int, default = 1, help = "Number of threads (default = 1).")
	parser.add_argument("-i", 
help = "Path to input sample (If omitted, the normal vcfs will be compared to one another).")
	parser.add_argument("-m", help = "Path to manifest of normals files (one file per line).")
	parser.add_argument("-o", help = "Path to output directory.")
	args = parser.parse_args()
	args, norm = checkArgs(args)
	normals, a, b = getNormals(args.m, args.o, args.allsamples)
	if norm == False and args.allsamples == False:
		vcfs, log = getSamplePairs(args.o, normals, args.i)
	elif args.allsamples == False:
		print("\tGetting all pairs of normal samples...")
		vcfs, log = getSamplePairs(args.o, normals)
	else:
		print("\tGetting all tumor:normal sample pairs...")
		vcfs, log = allSamplePairs(args.o, normals, a, b)
	func = partial(identifySample, norm, args.o, log)
	l = len(vcfs)
	pool = Pool(processes = args.t)
	print(("\tComparing vcf to normals with {} threads...\n").format(args.t))
	for x in pool.imap_unordered(func, vcfs):
		l -= 1
		if x[0] == False:
			print(("\t[Warning] Comparison between {} and {} failed.").format(x[1], x[2]))
		else:		
			print(("\tComparison between {} and {} successful. {} samples remaining.").format(x[1], x[2], l))
	print(("\tFinished. Runtime: {}\n").format(datetime.now()-start))

if __name__ == "__main__":
	main()
