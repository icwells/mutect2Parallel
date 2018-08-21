'''This script contains functions for adding readgroups to bam files, as well as indexing and extracting readgroups'''

import os
from unixpath import *

class Sample():
	# Stores data for managing sample progress
	def __init__(self):
		self.Name = ""
		self.ID = ""
		self.Step = ""
		self.Status = ""
		self.Status = ""
		self.Output = ""
		self.Private = None
		self.Bed = None
		self.Bam = None
		self.Input = None
		self.Unfiltered = None

	def update(self, sample, name, step, status, outfile):
		# Sorts and updates entry with additional status update
		if not self.Name:
			self.Name = sample
		if not self.ID:
			self.ID = name
		if step == "comparison":
			self.Step = step
			self.Output = outfile
			self.Status = status
		elif step == "filtering_covB" and self.Step != "comparison":
			self.Step = step
			self.Output = outfile
			self.Status = status
		elif step == "filtering_germline" and self.Step == "mutect":
			self.Step = step
			self.Output = outfile
			self.Status = status
		elif step == "mutect":
			if status == "complete" or self.Status == "starting":
				self.Unfiltered = outfile
				self.Step = step
				self.Output = outfile
				self.Status = status
		if getExt(outfile) == "bam":
			self.Bam = outfile

	def updateStatus(self, status, step = None, outfile = None, unfilt = False):
		# Updates current status of sample
		self.status = status
		if step:
			self.Step = step
		if outfile:
			self.Output = outfile
			if getExt(outfile) == "bam":
				self.Bam = outfile
		if unfilt == True:
			self.Unfiltered = outfile

#-------------------------------commonfunctions----------------------------------------

def getDelim(line):
	# Returns delimiter
	for i in ["\t", ",", " "]:
		if i in line:
			return i
	print("\n\t[Error] Cannot determine delimeter. Check file formatting. Exiting.\n")
	quit()

def getStatus(log):
	# Returns true if tool returns success
	status = False
	ret = False
	with open(log, "r") as l:
		for line in l:
			if "Tool returned:" in line:
				status = True
			elif status == True:
				# Get exit status
				if "SUCCESS" in line:
					ret = True
				break
	return ret

def appendLog(conf, s):
	# Appends checkpoint status to log file
	if s.Step == "mutect" and self.Status == "starting":
		# Record infile instead of outfile
		out = s.Infile
	elif "isec1" in s.Step and s.Private:
		# Record private vcf
		out = s.Private
	else:
		out = self.Output
	with open(conf["log"], "a") as l:
			l.write(("{}\t{}\t{}\t{}\t{}\n").format(s.Sample,s.ID, s.Step, s.Status, out))

def getOpt(conf, cmd):
	# Adds common flags to command
	if "bed" in conf.keys():
		cmd += (" -L {}").format(conf["bed"])
	if "germline" in conf.keys():
		cmd += (" --germline-resource {} --af-of-alleles-not-in-resource {}").format(
													conf["germline"], conf["af"])
	return cmd

def getSample(fname):
	# Returns sample name (ie raw filename)
	s = os.path.split(fname)[1]
	if "-" in s:
		# Remove group ID
		return s[s.find("-")+1:s.find(".")]
	else:
		return s[:s.find(".")]

def checkOutput(outdir, normal, prnt = True):
	# Checks for output log file and reads if present
	first = True
	done = {}
	log = outdir + "mutectLog.txt"
	if prnt == True:
		print("\tChecking for previous output...")
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	if os.path.isfile(log):
		with open(log, "r") as f:
			for line in f:
				if first == False and line.strip():
					line = line.strip().split("\t")
					if len(line) == 5:
						if line[0] not in done.keys():
							# Initialize new sample entry
							done[line[0]] = Sample()
						done[line[0]].update(line[0], line[1], line[2], line[3], line[4])
				else:
					# Skip header
					first = False
	else:
		with open(log, "w") as f:
			# Initialize log file and record normal file
			f.write("Sample\tName\tStep\tStatus\tOutput\n")
			f.write("N\t{}\tnormal\tcomplete\t{}".format(getFileName(normal), normal))
	return log, done

def configEntry(conf, arg, key):
	# Returns dict with updated arg entry
	if not arg:
		print(("\n\t[Error] Please specify {}. Exiting.\n").format(arg))
		quit()
	else:
		conf[key] = arg
	return conf

def getOptions(conf, line):
	# Returns config dict with updated options
	val = None
	s = line.split("=")
	if len(s) == 2:
		target = s[0].strip()
		val = s[1].strip()
	if val:
		if target == "reference_genome":
			conf["ref"] = val
		elif target == "bed_annotation":
			conf["bed"] = val
			checkFile(conf["bed"])
		elif target == "output_directory":
			if val[-1] != "/":
				val += "/"
			conf["outpath"] = val
		elif target == "GATK_jar":
			conf["gatk"] = val
			checkFile(conf["gatk"])
		elif target == "Picard_jar":
			conf["picard"] = val
			checkFile(conf["picard"])	
		elif target == "SnpSift_jar":
			conf["snpsift"] = val
			checkFile(conf["snpsift"])	
		elif target == "normal_panel":
			conf["pon"] = val
			checkFile(conf["pon"])
		elif target == "germline_resource":
			conf["germline"] = val
			checkFile(conf["germline"])
		elif target == "allele_frequency":
			conf["af"] = val
		elif target == "mutect_options":
			# Extract directly from line in case options have an equals sign
			conf["mo"] = line[line.find("=")+1:]
			conf["mo"] = conf["mo"].strip()		
		elif target == "contaminant_estimate":
			# for filterVCFs only
			conf["contaminant"] = val
			checkFile(conf["contaminant"])
		# Get filtereing parameters
		elif target == "qual":
			conf["qual"] = int(val)
		elif target == "min_covA":
			conf["min_covA"] = int(val)
		elif target == "min_reads_strand":
			conf["min_reads_strand"] = int(val)
		elif target == "min_reads_alt":
			conf["min_reads_alt"] = int(val)
		elif target == "min_covB":
			conf["min_covB"] = int(val)
		elif target == "max_altB":
			conf["max_altB"] = int(val)
		elif target == "max_prop_altB":
			conf["max_prop_altB"] = float(val)
		elif target == "max_covN":
			conf["max_covN"] = int(val)
		elif target == "min_freq_altN":
			conf["min_freq_altN"] = float(val)
		elif target == "min_reads_altN":
			conf["min_reads_altN"] = int(val)
	return conf

def getConf(infile):
	# Stores runtime options
	batch = []
	opt = True
	store = False
	conf = {"ref":None}
	print("\n\tReading config file...")
	with open(infile, "r") as f:
		for line in f:
			if "Sample Batch Script" in line:
				opt = False
			elif opt == True and line.strip():
				# Store options
					conf = getOptions(conf, line)
			else:
				if "#!/" in line:
					store = True
				if store == True:
					# Store batch script
					batch.append(line)
	# Check for critical errors
	if not os.path.isfile(conf["ref"]):
		print("\n\t[Error] Genome fasta not found. Exiting.\n")
		quit()
	if not os.path.isdir(conf["outpath"]):
		os.mkdir(conf["outpath"])
	if "germline" in conf.keys() and "af" not in conf.keys():
		print("\n\t[Error] Please supply an allele frequency when using a germline estimate. Exiting.\n")
		quit()	
	return conf, batch
