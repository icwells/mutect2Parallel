'''This script will run a pair of input tumor bam files through mutect2 in parallel'''

import os
from argparse import ArgumentParser
from datetime import datetime
from functools import partial
from multiprocessing import Pool, cpu_count
from runMutect import *

class Sample():
	# Stores data for managing sample progress
	def __init__(self, name, mostRecent, outfile=None):
		self.ID = name
		self.Status = "starting"
		if "failed" not in mostRecent:
			self.Status = mostRecent
		self.Output = outfile
		self.Bam = None
		self.Input = None

	def __update__(self, name, mostRecent, outfile):
		# Sorts and updates entry with additional status update
		if self.Status == "completed":
			pass
		elif self.Status == "":
			if "failed" not in mostRecent:
				self.Status = mostRecent
		elif mostRecent == "completed":
			self.Status = mostRecent
			self.Output = outfile
		elif mostRecent == "filtered" and self.Status == "mutect":
			self.Status = mostRecent
			self.Output = outfile
		elif mostRecent == "mutect" and self.Status == "starting":
			self.Status = mostRecent
			self.Output = outfile

#-----------------------------------------------------------------------------

def checkOutput(outdir):
	# Checks for output log file and reads if present
	first = True
	done = {}
	log = outdir + "mutectLog.txt"
	print("\tChecking for previous output...")
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	if os.path.isfile(log):
		with open(log, "r") as f:
			for line in f:
				if first == False and line.strip():
					line = line.strip().split("\t")
					if line[1] == "completed":
						done[line[0]] = Sample(line[0], line[1])
					if len(line) == 3:
						if line[0] in done.keys():
							done[line[0]].__update__(line[0], line[1], line[2])
						else:
							done[line[0]] = Sample(line[0], line[1], line[2])
				else:
					# Skip header
					first = False
	else:
		with open(log, "w") as f:
			# Initialize log file
			f.write("Filename\tStatus\tOutput\n")
	return log, done

def configEntry(conf, arg, key):
	# Returns dict with updated arg entry
	if not arg:
		print(("\n\t[Error] Please specify {}. Exiting.\n").format(arg))
		quit()
	else:
		conf[key] = arg
	return conf	

def getConfig(args):
	# Returns arguments as dict
	conf = {}
	conf["bamout"] = args.bamout
	if args.o[-1] != "/":
		args.o += "/"
	conf = configEntry(conf, args.s, "sample")
	conf = configEntry(conf, args.x, "tumor1")
	conf = configEntry(conf, args.y, "tumor2")
	conf = configEntry(conf, args.c, "normal")
	conf = configEntry(conf, args.r, "reference")
	conf = configEntry(conf, args.o, "outpath")
	if args.bed:
		conf["bed"] = args.bed
	if args.gatk:
		conf["gatk"] = args.gatk
	if args.picard:
		conf["picard"] = args.picard
	if args.p:
		conf["pon"] = args.p
	if args.g:
		if not args.af:
			print("\n\t[Error] Please supply an allele frequency when using a germline estimate. Exiting.\n")
			quit()
		else:
			conf["germline"] = args.g
			conf["af"] = args.af
	if args.e:
		conf["contaminant"] = args.e	
	return conf

def main():
	starttime = datetime.now()
	parser = ArgumentParser(description = "This script will call MuTect2 on a given \
list of input files. Be sure that pysam is installed and that bcftools is in your PATH.")
	parser.add_argument("--bamout", action = "store_true", default = False,
help = "Indicates that mutect should also generate bam output files.")
	parser.add_argument("-s", help = "Sample name (required).")
	parser.add_argument("-x", help = "Path to first tumor bam (required).")
	parser.add_argument("-y", help = "Path to second tumor bam (required).")
	parser.add_argument("-c", help = "Path to normal/control bam (required).")
	parser.add_argument("-r", help = "Path to reference genome (required).")
	parser.add_argument("-o", help = "Path to output directory (required).")
	parser.add_argument("--bed", help = "Path to bed annotation.")
	parser.add_argument("--gatk", help = "Path to gatk jar (if using).")
	parser.add_argument("--picard", help = "Path to picard jar (if using).")
	parser.add_argument("-p", help = "Path to panel of normals.")
	parser.add_argument("-g", help = "Path to germline resource.")
	parser.add_argument("--af", help = "Estimated allele frequency (required if using a germline resource).")
	parser.add_argument("-e", help = "Path to contmination estimate vcf.")
	args = parser.parse_args()
	conf = getConfig(args)
	log, samples = checkOutput(conf["outpath"])
	conf["log"] = log
	pool = Pool(processes = 2)
	func = partial(submitFiles, conf, samples)
	filtered = []
	# Call mutect
	print(("\n\tCalling mutect2 on {}....").format(conf["sample"]))
	for x in pool.imap_unordered(func, [conf["tumor1"], conf["tumor2"]]):
		if "failed" in x.Status:
			print(("\n\tFailed to run {}.").format(x.ID))
		else:		
			print(("\n\t{} has finished filtering.").format(x.ID))
			filtered.append(x.Output)
	pool.close()
	pool.join()
	if len(filtered) == 2:
		# Compare output
		print(("\n\tComparing filterd VCFs from {}...").format(conf["sample"]))
		status = compareVCFs(conf, filtered)
		if status == True:
			# Record finished samples
			outfile = conf["outpath"] + conf["sample"] + ".csv"
			with open(conf["log"], "a") as l:
				l.write(("{}\tcompleted\t{}\n").format(conf["sample"], outfile))
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
