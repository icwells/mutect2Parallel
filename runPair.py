'''This script will run a pair of input tumor bam files through mutect2 in parallel'''

import os
from argparse import ArgumentParser
from datetime import datetime
from functools import partial
from multiprocessing import Pool, cpu_count
from subprocess import Popen
from shlex import split
from bamUtil import *
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

#-------------------------------Mutect----------------------------------------

def getStatus(log):
	# Returns true if tool returns success
	status = False
	with open(log, "r") as l:
		for line in l:
			if "Tool returned:" in line:
				status = True
			elif status == True:
				# Get exit status
				if "SUCCESS" in line:
					return True
				else:
					return False

def callMutect(cmd, name, outfile):
	# Calls Mutect with given command
	print(("\tCalling mutect on {}...").format(name))
	# Make log file
	log = outfile.replace("vcf", "stdout")
	with open(log, "w") as dn:
		try:
			mt = Popen(split(cmd), stdout = dn, stderr = dn)	
			mt.communicate()
		except:
			print(("\t[Error] Could not call MuTect2 on {}").format(name))
			return None
	if getStatus(log) == True:
		print(("\t{} has completed mutect analysis.").format(name))
		return outfile
	else:
		return None

def getOpt(conf, cmd):
	# Adds common flags to command
	if "bed" in conf.keys():
		cmd += (" -L {}").format(conf["bed"])
	if "germline" in conf.keys():
		cmd += (" --germline-resource {} --af-of-alleles-not-in-resource {}").format(
													conf["germline"], conf["af"])
	return cmd

def submitSample(infile, conf, s, name):
	# Builds mutect command
	if "picard" in conf.keys():
		_, control = checkRG(conf["normal"], s.ID, conf["picard"])
		tumorname, bam = checkRG(infile, s.ID, conf["picard"])
	else:
		_, control = checkRG(conf["normal"], s.ID)
		tumorname, bam = checkRG(infile, name)
	if not control or not tumorname or not bam:
		s.Status = "failed-addingReadGroups"
		appendLog(conf, s)
		return s
	# Assemble command
	if "gatk" in conf.keys():
		# Format command for calling gatk jar
		cmd = ("java -jar {} Mutect2 -R {} ").format(conf["gatk"], conf["reference"])
	else:
		# Format command for colling gatk from path
		cmd = ("gatk Mutect2 -R {} ").format(conf["reference"])
	cmd += ("--tumor-sample {} -I {} -I {} --output {}").format(tumorname, 
										bam, conf["normal"], s.Output)
	if "bamout" in conf.keys():
		s.Bam = s.Output[:s.Output.find(".")] + ".Mutect2.bam"
		cmd += (" --bamout {}").format(s.Bam)
	if "pon" in conf.keys():
		cmd += (" --panel-of-normals {}").format(conf["pon"])
	if "mo" in conf.keys():
		cmd += " " + conf["mo"]
	cmd = getOpt(conf, cmd)
	# Call mutect for control and tumor
	res = callMutect(cmd, name, s.Output)
	if res:
		# Record finished sample
		s.Output = res
		s.Status = "mutect"
	else:
		s.Status = "failed-mutect"
	appendLog(conf, s)
	return s

#-----------------------------------------------------------------------------

def appendLog(conf, s):
	# Appends checkpoint status to log file
	with open(conf["log"], "a") as l:
			l.write(("{}\t{}\t{}\n").format(s.ID, s.Status, s.Output))

def getSample(fname):
	# Returns sample name (ie raw filename)
	s = os.path.split(fname)[1]
	if "-" in s:
		# Remove group ID
		return s[s.find("-")+1:s.find(".")]
	else:
		return s[:s.find(".")]

def submitFiles(conf, samples, infile):
	# Calls MuTect2 serially over input files
	name = getSample(infile)
	# Get sample info
	if name in samples.keys():
		s = samples[name]
	else:
		s = Sample(name, "starting", conf["outpath"] + name + ".vcf")
	s.Input = infile
	if s.Status == "starting":
		s = submitSample(infile, conf, s, name)
	return s

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
	if args.nofilter == True and args.filter == True:
		print("\n\t[Error] Please specify only one of filter/nofilter. Exiting.\n")
		quit()
	conf["bamout"] = args.bamout
	conf["filter"] = args.filter
	conf["nofilter"] = args.nofilter
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
	if args.mo:
		conf["fmo"] = args.mo
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
	parser.add_argument("--mo", help = "Additional mutect options in quotes")
	args = parser.parse_args()
	conf = getConfig(args)
	log, samples = checkOutput(conf["outpath"])
	conf["log"] = log
	if conf["filter"] == False:
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
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
