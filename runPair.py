'''This script will run a pair of input tumor bam files through mutect2 in parallel'''

import os
from argparse import ArgumentParser
from datetime import datetime
from functools import partial
from multiprocessing import Pool, cpu_count
from datetime import datetime
from subprocess import Popen
from shlex import split
from bamUtil import *

class Sample():
	# Stores data for managing sample progress
	def __init__(name, mostRecent, outfile):
		self.ID = name
		self.Status = mostRecent
		self.Output = outfile	

	def __update__(self, name, mostRecent, outfile):
		# Sorts and updates entry with additional status update
		if self.Status == "completed":
			pass
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

def getTotal(vcf):
	# Returns total number of content lines from vcf
	count = 0
	if os.path.isfile(vcf):
		with open(vcf, "r") as f:
			for line in f:
				if line[0] != "#":
					count += 1
	return count

def intersect(outpath, cmd, vcfs):
	# Calls bcftools to get intersecting rows and summarizes output
	cmd = cmd.format(outpath)
	try:
		bcf = Popen(split(cmd))
		bcf.communicate()
	except:
		print(("\t[Error] Could not call bcftools with {}").format(cmd))
		return 0
	# Number of unique variants to each sample and number of shared
	a = getTotal(outpath + "/0000.vcf")
	b = getTotal(outpath + "/0001.vcf")
	c = getTotal(outpath + "/0002.vcf")
	# Get percentage of similarity
	try:
		sim = c/(a+b+c)
	except ZeroDivisionError:
		sim = 0.0
	sa = vcfs[0][vcfs[0].rfind("/")+1:vcfs[0].find(".")]
	sb = vcfs[1][vcfs[1].rfind("/")+1:vcfs[1].find(".")]
	with open(outpath + ".csv", "w") as output:
		output.write("SampleA,SampleB,#PrivateA,#PrivateB,#Common,Similarity\n")
		output.write(("{},{},{},{},{},{}\n").format(sa, sb, a, b, c, sim))
	return 1

def compareVCFs(conf, vcfs):
	# Calls gatk and pyvcf to filter and compare mutect output
	ret = False
	done = 0
	outpath = conf["outpath"] + conf["sample"]
	# Call bftools on all results
	cmd = ("bcftools isec {} {}").format(vcfs[0], vcfs[1])
	cmd += " -p {}"
	done += intersect(outpath, cmd, vcfs)
	# Call bftools on passes
	outpath += "_PASS"
	done += intersect(outpath, cmd + " -f .,PASS", filtered)
	if done == 2:
		ret = True
	return ret

#-------------------------------Mutect----------------------------------------

def filterCalls(cmd, vcf):
	# Calls gatk to filter mutect calls
	outfile = vcf[:vcf.find(".")] + ".filtered.vcf"
	log = outfile.replace("vcf", "stdout")
	cmd += ("-V {} -O {}").format(vcf, outfile)
	with open(log, "w") as l:
		try:
			fmc = Popen(split(cmd), stdout = l, stderr = l)
			fmc.communicate()
		except:
			print(("\t[Error] Could not call FilterMutectCalls on {}").format(vcf))
			return ""
	return bgzip(outfile)

def callMutect(cmd, name):
	# Calls Mutect with given command
	print(("\tCalling mutect on {}...").format(name))
	# Make log file
	log = outfile.replace("vcf", "stdout")
	with open(log, "w") as dn:
		try:
			mt = Popen(split(cmd), stdout = dn, stderr = dn)	
			mt.communicate()
		except:
			print(("\t[Error] Could not call MuTect on {}").format(name))
			return ""
	with open(log, "r") as dn:
		status = False
		# Make sure mutect completed successfully
		for line in dn:
			if "Tool returned:" in line:
				status = True
			elif status == True:
				# Get exit status
				if "SUCCESS" in line:
					print(("\t{} has completed mutect analysis.").format(name))
					return outfile
				else:
					return ""

def submitFiles(conf, samples, infile):
	# Calls MuTect2 serially over input files
	name = infile[infile.rfind("/")+1:infile.find(".")]
	if name in samples.keys():
		s = sample[name]
	else:
		s = Sample(name, "starting", conf["outpath"] + name + ".vcf")
	if s.Status == "starting":
		if "picard" in conf.keys():
			_, control = checkRG(conf["normal"], s.ID, conf["picard"])
			tumorname, bam = checkRG(infile, name, conf["picard"])
		else:
			_, control = checkRG(conf["normal"], s.ID)
			tumorname, bam = checkRG(infile, name)
		# Assemble command
		if "gatk" in conf.keys():
			# Format command for calling gatk jar
			cmd = ("java -jar {} Mutect2 -R {} ").format(conf["gatk"], conf["ref"])
			filt = ("java -jar {} FilterMutectCalls ").format(conf["gatk"])
		else:
			# Format command for colling gatk from path
			cmd = ("gatk Mutect2 -R {} ").format(conf["ref"])
			filt = "gatk FilterMutectCalls "
		cmd += ("--tumor-sample {} -I {} -I {} --output {}").format(tumorname, bam, conf["normal"], s.Output)

		############## Add bed file ################################

		# Call mutect for control and tumor
		res = callMutect(cmd, name)
		if res:
			# Record finished sample
			s.Output = res
			s.Status = "mutect"
			with open(conf["log"], "a") as l:
				l.write(("{}\t{}\t{}\n").format(s.ID, s.Status, s.Output))
		else:
			s.Status = "failed-mutect"
	if s.Status == "mutect":
		filtered = filterCalls(filt, s.Output)
		if filtered:
			# Record filtered reads
			s.Output = filtered
			s.Status = "filtered"
			with open(conf["log"], "a") as l:
				l.write(("{}\t{}\t{}\n").format(s.ID, s.Status, s.Output))
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
					if len(line) == 3:
						if line[0] in done.keys():
							done[line[0]].__udpate__(line)
						else:
							done[line[0]] = Sample(line)
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
	if args.regions:
		conf["regions"] = args.regions
	return picard

def main():
	starttime = datetime.now()
	parser = ArgumentParser(description = "This script will call MuTect2 on a given \
list of input files. Be sure that pysam is installed and that bcftools is in your PATH.")
	parser.add_argument("-s", help = "Sample name (required).")
	parser.add_argument("-x", help = "Path to first tumor bam (required).")
	parser.add_argument("-y", help = "Path to second tumor bam (required).")
	parser.add_argument("-c", help = "Path to normal/control bam (required).")
	parser.add_argument("-r", help = "Path to reference genome (required).")
	parser.add_argument("-o", help = "Path to output directory (required).")
	parser.add_argument("-bed", help = "Path to bed annotation.")
	parser.add_argument("-gatk", help = "Path to gatk jar (if using).")
	parser.add_argument("-picard", help = "Path to picard jar (if using).")
	parser.add_argument("-a", help = "Path to active regions file.")
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
			print(("\n\t{} has completed.").format(x.ID))
			filtered.append(x.Output)
	pool.close()
	pool.join()
	if len(filtered) == 2:
		# Compare output
		print(("\n\tComparing filterd VCFs from {}...").format(conf["sample"]))
		status = compareVCFs(conf, filtered)
		if status == True:
			# Record finished samples
			with open(conf["log"], "a") as l:
				l.write(("{}\tcompleted\n").format(conf["sample"]))
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
