'''Filters previous output to reflect only current target variants.'''

from argparse import ArgumentParser
from datetime import datetime
import os
import pandas as pd
import unixpath

class Variant():

	def __init__(self, line):
		self.patient = line["Patient"]
		self.shared	= line["Shared"]
		self.chr = line["Chr"]
		self.start = line["Start"]
		self.end = line["End"]
		self.ref = line["REF"]
		self.alt = line["ALT"]

	def equals(self, line):
		# Returns true if line cantains same variant
		ret = True
		if self.patient != line["Patient"]:
			ret = False
		elif self.shared != line["Shared"]:
			ret = False
		elif self.chr != line["Chr"]:
			ret = False
		elif self.start != line["Start"]:
			ret = False
		elif self.end != line["End"]:
			ret = False
		elif self.ref != line["REF"]:
			ret = False
		elif self.alt != line["ALT"]:
			ret = False
		return ret

class VariantsFilter():

	def __init__(self, args):
		unixpath.checkFile(args.i)
		unixpath.checkFile(args.v)
		self.infile = args.i
		self.outfile = ""
		self.variants = []
		self.xlsx = args.v
		self.__readVariants__(args)

	def __readVariants__(self, args):
		# Loads variants from xlsx file
		print("\n\tReading confirmed variants...")
		outdir = os.path.split(self.infile)[0]
		if args.ampliseq:
			sheet = "Table 3S_SNVs validation"
			name = "ampliseqVariants.csv"
		elif args.relaxed:
			sheet = "Table 4S_S_ list varinats (SNVs"
			name = "relaxedVariants.csv"
		elif args.stringent:
			sheet = "Table 5S_R_ list variants(SNVs)"
			name = "stringentVarinats.csv"
		# Store outfile name
		self.outfile = os.path.join(outdir, name)
		wb = pd.read_excel(args.v, sheet)
		for _, i in wb.iterrows():
			self.variants.append(Variant(i))

	def filterVariants(self):
		# Filters variants in infile
		print("\tFiltering input variants...")
		df = pd.read_csv(self.infile, sep = ",", header = 0)
		for idx, i in df.iterrows():
			if len(i["REF"]) != 1 or i["REF"] == "-":
				df = df.drop(index = idx)
			else:
				match = False
				for v in self.variants:
					if v.equals(i):
						match = True
						break
				if not match:
					df = df.drop(index = idx)
		df.to_csv(self.outfile)
		print(("\tIdentified {} variants.").format(df.shape[0]))

def main():
	start = datetime.now()
	parser = ArgumentParser("Filters previous output to reflect only current target variants.")
	parser.add_argument("--ampliseq", action = "store_true", default = False, help = "Filter ampliseq variants.")
	parser.add_argument("-i", help = "Path to mutect2 variants csv.")
	parser.add_argument("--relaxed", action = "store_true", default = False, help = "Filter relaxed variants.")
	parser.add_argument("--stringent", action = "store_true", default = False, help = "Filter stringent variants.")
	parser.add_argument("-v", help = "Path to xlsx file of approved variants.")
	f = VariantsFilter(parser.parse_args())
	f.filterVariants()
	print(("\tTotal runtime: {}\n").format(datetime.now() - start))

if __name__ == "__main__":
	main()
