'''This script will run a pair of input tumor bam files through mutect2 in parallel'''

import os
from argparse import ArgumentParser
from datetime import datetime
from functools import partial
from multiprocessing import Pool, cpu_count
from subprocess import Popen
from shlex import split
from commonUtil import *

#-------------------------------Mutect----------------------------------------

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

def submitSample(infile, conf, s, name):
	# Builds mutect command
	if "picard" in conf.keys():
		_, control = checkRG(conf["normal"], s.ID, conf["picard"])
		tumorname, bam = checkRG(infile, s.ID, conf["picard"])
	else:
		_, control = checkRG(conf["normal"], s.ID)
		tumorname, bam = checkRG(infile, name)
	if not control or not tumorname or not bam:
		s.Step = "addingReadGroups"
		s.Status = "failed"
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
	s.Status = "mutect"
	if res:
		# Record finished sample
		s.Output = res
		s.Status = "passed"
	else:
		s.Status = "failed"
	appendLog(conf, s)
	return s

def submitFiles(conf, samples, infile):
	# Creates sample entry and calls MuTect2 over input files
	name = getSample(infile)
	# Get sample info
	if name in samples.keys():
		s = samples[name]
	else:
		s = Sample(name, "preflight", "starting", conf["outpath"] + name + ".vcf")
	s.Input = infile
	if s.Status == "starting":
		s = submitSample(infile, conf, s, name)
	return s

#-----------------------------------------------------------------------------

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
	if args.mo:
		conf["mo"] = args.mo
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
	parser.add_argument("--mo", help = "Additional mutect options in quotes (these will not be checked for errors).")
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
		if x.Status == "failed":
			print(("\n\tFailed to run {}.").format(x.ID))
		else:		
			print(("\n\t{} has finished mutect.").format(x.ID))
			filtered.append(x.Output)
	pool.close()
	pool.join()
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
