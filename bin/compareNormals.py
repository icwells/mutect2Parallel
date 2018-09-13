'''This script will call bcftools isec to identify any mis-labled samples'''

from datetime import datetime
from argparse import ArgumentParser
from itertools import combinations
from multiprocessing import Pool, cpu_count
from functools import partial
from unixpath import *
from mutect2Parallel import getManifest
from commonUtil import *

def identifySample(outdir, log, vcfs):
	# Calls bcftools isec for each set of normals vs vcf
	v = getParent(vcfs[0])
	n = getParent(vcfs[1])
	outpath = outdir + n
	a = bcfIsec(outpath, vcfs)
	if a is not None:
		b = getTotal(outpath + "/0001.vcf")
		c = getTotal(outpath + "/0002.vcf")
		try:
			sim = c/(a+b+c)
		except ZeroDivisionError:
			sim = 0.0
		with open(log, "a") as out:
			out.write(("{},{},{},{},{},{},{:.2%}\n").format(v, n, a, b, c, sim))
		return [True, v, n]
	else:
		return [False, v, n]

def getSomaticVariants(args):
	# Calls platypus to remove germline variants
	log = args.o + getFileName(args.i) + ".txt"
	outfile = args.o + getFileName(args.i) + ".vcf"
	cmd = ("python {}platypus callVariants --nCPU={} --bamFiles={} ").format(args.p, args.t, args.i)
	cmd += ("--refFile={} --output={} --logFileName={}").format(args.r, outfile, log)
	cmd += " --filterReadPairsWithSmallInserts=,0"
	print("\tFiltering bam file...")
	res = runProc(cmd)
	if res == False:
		printFatal("Bam file failed platypus filtering")
	return tabix(outfile)

def getSamplePairs(outdir, normals, vcf = None):
	# Returns pairs of samples to compare
	if vcf:
		vcfs = []
		log = outdir + getFileName(vcf) + "Comparison.csv"
		for i in normals:
			# Pair input vcf with each normal vcf
			vcfs.append([vcf, i])
	else:
		log = outdir + "normalsComparison.csv"
		vcfs = combinations(normals, 2)
	with open(log, "w") as out:
		# Initialize log file
		out.write("SampleA, SampleB, PrivateA, PrivateB, Common, %Similarity\n")
	return vcfs, log

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
	args.o = checkDir(args.o)
	args.p = checkDir(args.p)
	if args.t > cpu_count():
		args.t = cpu_count()
	return args, norm

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will call bcftools isec to compare input samples.")
	parser.add_argument("-t", type = int, default = 1, help = "Number of threads (default = 1).")
	parser.add_argument("-i", 
help = "Path to input sample (If omitted, the normal vcfs will be compared to one another).")
	parser.add_argument("-m", help = "Path to manifest of normals files.")
	parser.add_argument("-o", help = "Path to output directory.")
	parser.add_argument("-p", default = "/packages/6x/platypus/0.8.1/",
 help = "Path to platypus package directory (default = /packages/6x/platypus/0.8.1/).")
	parser.add_argument("-r", default = "/home/dmalload/storage/DCIS/temp_storage/GRCh37-lite.fa",
 help = "Path to reference genome (default = /home/dmalload/storage/DCIS/temp_storage/GRCh37-lite.fa).")
	args = parser.parse_args()
	args, norm = checkArgs(args)
	print()
	normals = getManifest(args.m, True)
	if norm == False:
		vcf = getSomaticVariants(args)
		vcfs, log = getSamplePairs(args.o, normals, vcf)
	else:
		print("\n\tGetting all pairs of normal samples...")
		vcfs, log = getSamplePairs(args.o, normals)
	func = partial(identifySample, args.o, log)
	l = len(vcfs)
	pool = Pool(processes = args.t)
	print(("\tComparing vcf to normals with {} threads...\n").format(args.t))
	for x in pool.imap_unordered(func, normals):
		l -= 1
		if x[0] == False:
			print(("\t[Warning] Comparison between {} and {} failed.").format(x[1], x[2]))
		else:		
			print(("\tComparison between {} and {} successful. {} samples remaining.").format(x[1], x[2], l))
	print(("\tFinished. Runtime: {}\n").format(datetime.now()-start))

if __name__ == "__main__":
	main()
