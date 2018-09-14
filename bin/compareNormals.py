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

def getSomaticVariants(args):
	# Calls platypus to remove germline variants
	if args.s:
		outfile = args.o + args.s + ".vcf"
		log = args.o + args.s + ".txt"
	else:
		outfile = args.o + getFileName(args.i) + ".vcf"
		log = args.o + getFileName(args.i) + ".txt"
	if not os.path.isfile(outfile) or os.path.getsize(outfile) < 2097152:
		# Proceed if output file does not exist or is it less than 2mb
		cmd = ("platypus callVariants --nCPU={} --bamFiles={} ").format(args.t, args.i)
		cmd += ("--refFile={} --output={} --logFileName={}").format(args.r, outfile, log)
		cmd += " --filterReadPairsWithSmallInserts=0"
		print("\tFiltering bam file...")
		print(cmd)
		res = runProc(cmd)
		if res == False:
			printFatal("Bam file failed platypus filtering")
	else:
		print("\tProceeding with existing vcf file.")
	return tabix(outfile)

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

def getNormals(infile, outdir):
	# Reads input manifest and copies files to outdir
	normals = []
	outpath = checkDir(outdir + "normals/", True)
	print("\n\tReading normals manifest...")
	with open(infile, "r") as f:
		for line in f:
			line = line.strip()
			if os.path.isfile(line):
				# Get outfile name with sample as part of file name
				outfile = outpath + getParent(line) + ".N.vcf"
				if os.path.isfile(outfile + ".gz"):
					# Check for existance of gzipped file
					outfile += outfile + ".gz"
				elif not os.path.isfile(outfile):
					# Copy file if it has not already been copied
					copy(line, outfile)
				normals.append(outfile)
			else:
				print(("\t[Warning] {} not found. Skipping.").format(line))
	return normals

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
		checkFile(args.i)
	else:
		norm = True
	checkFile(args.m)
	checkFile(args.r)
	args.o = checkDir(args.o, True)
	args.p = checkDir(args.p)
	if args.t > cpu_count():
		args.t = cpu_count()
	return args, norm

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will call bcftools isec to compare input samples. \
Make sure platypus is loaded in a module or in your PATH if supplying an input bam file.")
	parser.add_argument("-t", type = int, default = 1, help = "Number of threads (default = 1).")
	parser.add_argument("-i", 
help = "Path to input sample (If omitted, the normal vcfs will be compared to one another).")
	parser.add_argument("-s", help = "Optional sample name (if not in input file name).")
	parser.add_argument("-m", help = "Path to manifest of normals files (one file per line).")
	parser.add_argument("-o", help = "Path to output directory.")
	parser.add_argument("-r", default = "/home/dmalload/storage/DCIS/temp_storage/GRCh37-lite.fa",
 help = "Path to reference genome (default = /home/dmalload/storage/DCIS/temp_storage/GRCh37-lite.fa).")
	args = parser.parse_args()
	args, norm = checkArgs(args)
	print()
	normals = getNormals(args.m, args.o)
	if norm == False:
		vcf = getSomaticVariants(args)
		vcfs, log = getSamplePairs(args.o, normals, vcf)
	else:
		print("\tGetting all pairs of normal samples...")
		vcfs, log = getSamplePairs(args.o, normals)
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
