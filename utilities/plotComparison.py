'''This script will make a scatter plot of the percent of similar variants 
between platypus and mutect2'''

import os
from datetime import datetime
from argparse import ArgumentParser
from sys import stderr
from unixpath import *
from plotFunctions import *

def mergeSums(m, p):
	# Returns paired lists of values such that each equivalent value has the same index
	x = []
	y = []
	for i in m.keys():
		if i in p.keys():
			x.append(m[i])
			y.append(p[i])
	return [x, y]

def getOutputSummary(val, infile, percent = True):
	# Retuns dict of target values
	first = True
	sim = {}
	with open(infile, "r") as f:
		for line in f:
			if first == False:
				s = line.strip().split(",")
				if len(s) > c.Max:
					i = s[c.ID]
					idx = i.find("_")
					if idx == 2 or idx == 3:
						# Trim preceding number
						i = i[idx+1:]
					if val == "s":
						if percent == True:
							n = float(s[c.Similarity].strip("%"))/100
						else:
							n = float(s[c.Similarity])
					elif val == "a":
						n = int(s[c.A])
					elif val == "b":
						n = int(s[c.B])
					elif val == "c":
						n = int(s[c.Common])
					# {id: n}
					sim[i] = n
			else:
				c = Columns(line.split(","))
				first = False
	return sim

def getSummaries(val, mutect, platypus):
	# Returns list of paired values to plot
	print("\n\tReading input files...")
	m = getOutputSummary(val, mutect, True)
	p = getOutputSummary(val, platypus, False)
	return mergeSums(m, p)

def checkArgs(args):
	# Checks for errors and returns output file name
	if not args.m or not args.p:
		print("\n\t[Error] Please specify mutect and platypus summary files. Exiting.\n", file = stderr)
		quit()
	checkFile(args.m)
	checkFile(args.p)
	if args.o:
		parent = getParent(args.o)
		parent = checkDir(parent)
		if args.o[args.o.rfind("."):] != ".svg":
			# Add svg extension
			args.o += ".svg"
	else:
		# Write output to same directory
		args.o = args.m.replace(".csv", ".svg")
	args.v = args.v.lower()
	good = False
	for i in ["s", "a", "b", "c"]:
		if args.v == i:
			good = True
			break
	if good == False:
		print("\n\t[Error] Please enter one of s, a, b, c for -v. Exiting.\n", file = stderr)
		quit()
	return args

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will make an svg scatter plot of the \
percent of similar variants between platypus and mutect2.")
	parser.add_argument("-v", default = "s", help = "Code for values to plot. \
Values include Percent similarity (default): s, private A: a, private b: b, common: c")
	parser.add_argument("-m", help = "Path to mutect2 summary file.")
	parser.add_argument("-p", help = "Path to platypus summary file.")
	parser.add_argument("-o", 
help = "Path to output svg (will be written to same directory by default).")
	args = parser.parse_args()
	args = checkArgs(args)
	points = getSummaries(args.v, args.m, args.p)
	plotSimilarity(args.o, args.v, points)
	print(("\tFinished. Runtime: {}\n").format(datetime.now()-start))

if __name__ == "__main__":
	main()
