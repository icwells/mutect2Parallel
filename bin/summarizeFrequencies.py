'''This script will summarize the frequencies of variants from the mutect2 and platypus pipelines.'''

import os
import gzip
from datetime import datetime
from argparse import ArgumentParser
from glob import glob
from unixpath import *
from commonUtil import *
from compareNormals import fatalError
from pipelineComparison import platypusPaths

class Variant():
	def __init__(self, chrom, coor):
		# Store data for each variant
		self.chr = chrom
		self.coor = coor
		self.mutectA = 0
		self.mutectB = 0
		self.mutectC = 0
		self.platypusA = 0
		self.platypusB = 0
		self.platypusC = 0

	def __str__(self):
		# Returns string of single line
		ret = ("{},{},{},").format(self.name, self.chr, self.coor)
		ret += ("{},{},{},").format(self.platypusC, self.platypusA, self.platypusB)
		ret += ("{},{},{}").format(self.mutectC, self.mutectA, self.mutectB)
		return ret

class Variants():
	def __init__(self, sample):
		# Defines class for storing variants for samples from one tumor
		self.name = sample
		self.output = {}
		# [[variant, chrom, coor, freq]]

	def __str__(self):
		# Converts lists to writable string
		ret = ""
		for i in self.output.keys():
			for j in self.output[i].keys():
				ret += ("{},{}\n").format(self.name, i)
		return ret

	def __readFrequency__(self, info):
		# Get population alt allele frequency from info column
		f = "NA"
		for i in info:
			if "POP_AF" in i:
				f = i.split("=")[-1]
				break
		return f

	def __assignFrequency__(self, column, chrom, coor, freq):
		# Assigns frequency to proper column of proper dict
		if chrom not in self.output.keys():
			self.output[chrom] = {}
		if coor not in self.output[chrom].keys():
			self.output[chrom][coor] = Variant(chrom, coor)
		if "m" in column:
			if "a" in column:
				self.output[chrom][coor].mutectA = freq
			elif "b" in column:
				self.output[chrom][coor].mutectB = freq
			elif "c" in column:
				self.output[chrom][coor].mutectC = freq
		elif "p" in column:
			if "a" in column:
				self.output[chrom][coor].platypusA = freq
			elif "b" in column:
				self.output[chrom][coor].platypusB = freq
			elif "c" in column:
				self.output[chrom][coor].platypusC = freq

	def __readVCF__(self, vcf, column):
		# Reads uncompressed vcf
		with open(vcf, "r") as f:
			for line in f:
				if line[0] != "#":
					s = line.strip().split("\t")
					f = self.__readFrequency__(s[7].split(";"))
					self.__assignFrequency__(column, s[0], s[1], f)

	def __readGZ__(self, vcf):
		# Reads gzippped vcf
		tag = ("#").encode()
		with gzip.open(vcf, "rb") as f:
			for line in f:
				if tag not in line:
					s = str(line).split("\t")
					chrom, coor = s[0], s[1]
					info = s[7].split(";")
					v.append(self.__readFrequency__(info))
					print(chrom)

	def __readVariants__(self, vcf):
		# Assigns file to appropriate reading function and returns dict
		ret = []
		if os.path.isfile(vcf):
			if ".gz" not in vcf:
				ret = self.__readVCF__(vcf)
			else:
				ret = self.__readGZ__(vcf)
		elif os.path.isfile(vcf + ".gz"):
			ret = self.__readGZ__(vcf + ".gz")
		return ret

	def getMutectVariants(self, path):
		# Reads variant frequencies from filtered mutect output files
		self.__readVariants__(p + "A.NAB.vcf", "ma")
		self.__readVariants__(p + "B.NAB.vcf", "mb")
		self.__readVariants__(p + "A.NAB.vcf", "mc")

	def getPlatypusVariants(self, paths):
		# Reads variant frequencies from filtered platypus output files
		self.platypusA = self.__readVariants__(paths["a"])
		self.platypusB = self.__readVariants__(paths["b"])
		self.platypusC = self.__readVariants__(paths["c"])

#-----------------------------------------------------------------------------

class Frequencies():
	def __init__(self, mdir, pdir, out):
		# Defines class for reading and storing variant frequencies
		self.header = "Patient,Chromosome,Coordinates,PlatypusCommon,PlatypusPrivateA,PlatypusPrivateB,Mutect2Common,Mutect2PrivateA,Mutect2PrivateB\n"
		self.mutectDir = mdir
		self.platypusDir = pdir
		self.outfile = out
		self.frequencies = {}

	def writeFrequencies(self):
		# Writes dict to file
		print("\tWriting frequencies to file...")
		with open(self.outfile, "w") as out:
			out.write(self.header)
			for i in self.frequencies.keys():
				out.write(str(self.frequencies[i]))

	def mutectFrequencies(self):
		# Stores frequencies from mutect2 output
		paths = glob(self.MutectDir + "*")
		print("\n\tGetting mutect2 variant frequencies...")
		for p in paths:
			if os.path.isdir(p):
				if p != "comparison" and p != "raw" and p != "unfiltered":
					full = getParent(p)
					# Drop leading number from name
					sample = full[full.find("_")+1:]
					v = Variants(sample)
					self.Frequencies[sample] = v.getMutectVariants(p)

	def platypusFrequencies(self):
		# Stores frequencies from mutect2 output
		paths = glob(self.MutectDir + "*")
		print("\tGetting platypus variant frequencies...")
		for p in paths:
			if os.path.isdir(p):
				paths = platypusPaths(p, "nab")
				full = getParent(p)
				if sample not in self.frequencies.keys():
					v = Variants(sample)
				self.Frequencies[sample] = v.getPlatypusVariants(paths)


def checkArgs(args):
	# Checks argument validity
	if not args.m or not args.p:
		fatalError("Please provide paths to mutect2 (-m) and platypus (-p) output directories")
	if not args.o:
		fatalError("Please provide name of output csv file")
	args.m = checkDir(args.m)
	args.p = checkDir(args.p)
	# Make sure parent directory exists
	_ = checkDir(os.path.split(args.o)[0])
	return args

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will summarize the frequencies of \
variants from the mutect2 and platypus pipelines.")
	parser.add_argument("-m", help = "Path to mutect output directory.")
	parser.add_argument("-p", help = "Path to platypus output directory.")
	parser.add_argument("-o", help = "Path to output csv.")
	args = parser.parse_args()
	args = checkArgs(args)
	f = Frequencies(args.m, args.p, args.o)
	f.mutectFrequencies()
	f.platypusFrequencies
	f.writeFrequencies()
	print(("\tFinished. Runtime: {}\n").format(datetime.now()-start))

if __name__ == "__main__":
	main()
