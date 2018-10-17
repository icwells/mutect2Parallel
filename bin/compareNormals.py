'''This script will call bcftools isec to identify any mis-labled samples'''

import os
from datetime import datetime
from argparse import ArgumentParser
from sys import stderr
from shutil import copy
from itertools import combinations
from multiprocessing import Pool, cpu_count
from unixpath import *
from commonUtil import *

class Finished():

	def __init__(self, infile):
		self.Log = infile
		self.A = {}
		self.B = {}
		self.N = {}
		self.__getFinished__()

	def __getFinished__(self):
		# Returns list of pairs which finished comparison
		first = True
		if not os.path.isfile(self.Log):
			print("\tMaking new log file...")
			with open(self.Log, "w") as out:
				# Initialize log file
				out.write("SampleType,Sample,Normal,PrivateSample,PrivateNormal,Common,%Similarity\n")
		else:
			with open(self.Log, "r") as f:
				for line in f:
					if first == False:
						s =  line.split(",")
						# Add to appropriate dict
						if s[0] == "A":
							if s[1] not in self.A.keys():
								self.A[s[1]] = []
							self.A[s[1]].append(s[2])
						elif s[0] == "B":
							if s[1] not in self.B.keys():
								self.B[s[1]] = []
							self.B[s[1]].append(s[2])
						elif s[0] == "N":
							if s[1] not in self.N.keys():
								self.N[s[1]] = []
							self.N[s[1]].append(s[2])
					else:
						first = False
			l = len(self.A) + len(self.B) + len(self.N)
			print(("\tIdentified {:,d} completed samples.").format(l))

	def inFinished(self, typ, vcf1, vcf2):
		# Reuturns true if vcfs pair is in done
		v = getFileName(vcf1)
		n = getFileName(vcf2)
		if typ == "A":
			if v in self.A.keys():
				if n in self.A[v]:
					return True
		elif typ == "B":
			if v in self.B.keys():
				if n in self.B[v]:
					return True
		elif typ == "N":
			if v in self.N.keys():
				if n in self.N[v]:
					return True
		return False

#-----------------------------------------------------------------------------

class VCFcomparison():

	def __init__(self, t, v, n, log, outdir, norm = None, sample = None):
		self.type = t
		self.vcf = v
		self.normal = n
		self.log = log
		self.v = getFileName(self.vcf)
		self.n = getFileName(self.normal)
		self.outdir = ""
		self.sample = sample
		self.__checkInput__()
		self.__getOutdir__(outdir, norm)

	def __getOutdir__(self, outdir, norm):
		# Gets outdir for bcftools
		if norm is not None:
			# Get outdir for platypus comparisons
			root = "isecSample"
			if norm == True:
				root = "isecNormals"
			self.outdir = ("{}{}/{}_{}").format(outdir, root, self.v, self.n)	
		else:
			# Get pipeline output
			self.outdir = outdir + self.sample + "/"

	def __checkInput__(self):
		# Makes sure input files exist
		self.vcf = checkGZ(self.vcf)
		self.normal = checkGZ(self.normal)

#-----------------------------------------------------------------------------

def getType(vcf):
	# Returns sample type from filename
	if ".A." in vcf:
		return "A"
	elif ".B." in vcf:
		return "B"
	else:
		return "N"

def compareSamples(v):
	# Calls bcftools isec for each set of normals vs vcf
	a = bcfIsec(v.outdir, [v.vcf, v.normal])
	if a is not None:
		b = getTotal(v.outdir + "/0001.vcf")
		c = getTotal(v.outdir + "/0002.vcf")
		try:
			sim = c/(a+b+c)
		except ZeroDivisionError:
			sim = 0.0
		with open(v.log, "a") as out:
			t = v.type
			vcf = v.v
			n = v.n
			if v.sample is not None:
				t = v.sample + "," + v.type
				vcf = v.vcf
				n = v.normal
			out.write(("{},{},{},{},{},{},{:.2%}\n").format(t, vcf, n, a, b, c, sim))
		return [True, v.v, v.n]
	else:
		return [False, v.v, v.n]

def allSamplePairs(outdir, normals, a, b):
	# Returns all pairs for a:normal and b:normal
	vcfs = []
	log = outdir + "allSamplesComparison.csv"
	done = Finished(log)
	for t in [a, b]:
		for i in t:
			typ = getType(i)
			for j in normals:
				if done.inFinished(typ, i, j) == False:
					# Append each a/b to normal pair
					c = VCFcomparison(typ, i, j, log, outdir, True)
					vcfs.append(c)
	return vcfs

def getSamplePairs(outdir, normals, vcf = None):
	# Returns pairs of samples to compare
	vcfs = []
	if vcf:
		log = outdir + getFileName(vcf) + "Comparison.csv"
		done = Finished(log)
		typ = getType(vcf)
		for i in normals:
			if done.inFinished(typ, vcf, i) == False:
				# Pair input vcf with each normal vcf
				c = VCFcomparison(typ, vcf, i, log, outdir, False)
				vcfs.append(c)
	else:
		log = outdir + "normalsComparison.csv"
		done = Finished(log)
		c = combinations(normals, 2)
		for i in c:
			typ = getType(i[0])
			if done.inFinished(typ, i[0], i[1]) == False:
				c = VCFcomparison(typ, i[0], i[1], log, outdir, True)
				vcfs.append(c)
	return vcfs

def checkVCF(inpath, outpath, stem):
	# Returns name of existing output file and copies if necessary
	outfile = None
	infile = ("{}/{}").format(inpath, stem)
	if os.path.isfile(infile):
		outfile = ("{}{}.{}").format(outpath, getParent(inpath), stem)
		if os.path.isfile(outfile + ".gz"):
			# Check for existance of gzipped file
			outfile += ".gz"
		elif not os.path.isfile(outfile):
			# Copy file if it has not already been copied
			copy(infile, outfile)
	else:
		print(("\t[Warning] {} not found. Skipping.").format(line), file=stderr)
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
			path = os.path.split(line)[0]
			# Get outfile name with sample as part of file name
			outfile = checkVCF(path, outpath, "N.vcf")
			if outfile:
				normals.append(outfile)
			if allsamples == True:
				afile = checkVCF(path, outpath, "A.vcf")
				if afile:
					tumora.append(afile)
				bfile = checkVCF(path, outpath, "B.vcf")
				if bfile:
					tumorb.append(bfile)
	return normals, tumora, tumorb

def fatalError(msg):
	# Prints meassage and exits
	print(("\n\t[Error] {}. Exiting.\n").format(msg), file=stderr)
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
		print("\tGetting all sample:normal pairs...")
		vcfs = getSamplePairs(args.o, normals, args.i)
	elif args.allsamples == False:
		print("\tGetting all pairs of normal samples...")
		vcfs = getSamplePairs(args.o, normals)
	else:
		print("\tGetting all tumor:normal sample pairs...")
		vcfs = allSamplePairs(args.o, normals, a, b)
	print(("\t{:,d} file pairs found.").format(len(vcfs)))
	l = len(vcfs)
	pool = Pool(processes = args.t)
	print(("\tComparing vcf to normals with {} threads...\n").format(args.t))
	for x in pool.imap_unordered(compareSamples, vcfs):
		l -= 1
		if x[0] == False:
			print(("\t[Warning] Comparison between {} and {} failed.").format(x[1], x[2]), flush = True)
		else:		
			print(("\tComparison between {} and {} successful. {:,d} sets remaining.").format(x[1], x[2], l), flush = True)
	print(("\tFinished. Runtime: {}\n").format(datetime.now()-start))

if __name__ == "__main__":
	main()
