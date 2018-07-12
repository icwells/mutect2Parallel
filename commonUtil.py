'''This script contains functions for adding readgroups to bam files, as well as indexing and extracting readgroups'''

import os
import pysam
from subprocess import Popen
from shlex import split

class Sample():
	# Stores data for managing sample progress
	def __init__(self):
		self.ID = ""
		self.Step = ""
		self.Status = ""
		self.Status = ""
		self.Output = ""
		self.Bam = None
		self.Input = None
		self.Unfiltered = None

	def update(self, name, step, status, outfile):
		# Sorts and updates entry with additional status update
		if not self.ID:
			self.ID = name
		if step == "comparison":
			self.Step = step
			self.Output = outfile
			self.Status = status
		elif step == "filtering" and self.Step == "mutect":
			self.Step = step
			self.Output = outfile
			self.Status = status
		elif step == "mutect":
			self.Unfiltered = outfile
			self.Step = step
			self.Output = outfile
			self.Status = status

#-------------------------------commonfunctions----------------------------------------

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

def appendLog(conf, s):
	# Appends checkpoint status to log file
	with open(conf["log"], "a") as l:
			l.write(("{}\t{}\t{}\t{}\n").format(s.ID, s.Step, s.Status, s.Output))

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

def checkOutput(outdir, prnt = True):
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
					if len(line) == 4:
						if line[0] not in done.keys():
							# Initialize new sample entry
							done[line[0]] = Sample()
						done[line[0]].update(line[0], line[1], line[2], line[3])
				else:
					# Skip header
					first = False
	else:
		with open(log, "w") as f:
			# Initialize log file
			f.write("Filename\tStep\tStatus\tOutput\n")
	return log, done

def checkFile(infile, err):
	# Exits if infile is not found
	if not os.path.isfile(infile):
		print(("\n\t[Error] {} not found. Exiting.\n").format(err))
		quit()

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
			checkFile(conf["bed"], "Bed file")
		elif target == "output_directory":
			if val[-1] != "/":
				val += "/"
			conf["outpath"] = val
		elif target == "GATK_jar":
			conf["gatk"] = val
			checkFile(conf["gatk"], "GATK jar")
		elif target == "Picard_jar":
			conf["picard"] = val
			checkFile(conf["picard"], "Picard jar")	
		elif target == "normal_panel":
			conf["pon"] = val
			checkFile(conf["pon"], "Panel of normals file")
		elif target == "germline_resource":
			conf["germline"] = val
			checkFile(conf["germline"], "Germline resource file")
		elif target == "allele_frequency":
			conf["af"] = val
		elif target == "mutect_options":
			# Extract directly from line in case options have an equals sign
			conf["mo"] = line[line.find("=")+1:]
			conf["mo"] = conf["mo"].strip()	
		elif target == "mutect_options":
			# for filterVCFs only
			conf["fmo"] = line[line.find("=")+1:]
			conf["fmo"] = conf["fmo"].strip()	
		elif target == "contaminant_estimate":
			# for filterVCFs only
			conf["contaminant"] = val
			checkFile(conf["contaminant"], "Contamination estimate file")
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


#-------------------------------bamfunctions----------------------------------------

def bgzip(vcf):
	# bgzip compresses filtered vcf files
	if os.path.isfile(vcf + ".gz"):
		gz = vcf + ".gz"
	else:	
		gz = pysam.tabix_index(vcf, seq_col=0, start_col=1, end_col=1, force=True)
	return gz

def getFastaIndex(ref):
	# Creates fasta index for reference genome
	print("\tGenerating reference fasta index...")
	pysam.faidx(ref)

def samIndex(bam):
	# Call samtools to index sam/bam file
	if not os.path.isfile(bam + ".bai"):
		print(("\tGenerating sam index for {}...").format(bam))
		pysam.index(bam)

def getTumorName(bam):
	# Gets relevent header data and extracts tumor sample name from bam file
	try:
		header = pysam.view("-H", bam)
		rg = header[header.find("@RG"):header.find("@PG")]
		name = rg[rg.find("SM:"):]
		return name[name.find(":")+1:name.find("\t")]
	except pysam.utils.SamtoolsError:
		return None

def addRG(bam, sid, picard):
	# Adds readgroups to bam files
	outfile = bam[:bam.rfind(".")] + ".withRG.bam"
	if picard:
		cmd = ("java -jar {} AddOrReplaceReadGroups I={} O={} ").format(picard, bam, outfile)
	else:
		cmd = ("picard AddOrReplaceReadGroups I={} O={}").format(bam, outfile)
	cmd += (" RGLB=lib1 RGPL=illumina RGPU={}").format(sid)
	try:
		print(("\tAdding read groups to {}").format(bam))
		with open(os.devnull, "w") as dn:
			arg = Popen(split(cmd), stdout = dn, stderr = dn)
			arg.communicate()
	except:
		print(("\t[Error] Could not add read groups to {}\n").format(bam))
		return None, None
	if outfile:
		name = getTumorName(outfile)
		return outfile, name

def checkRG(bam, sid, picard=None):
	# Adds read groups, creates bam index, and returns read group name
	idx = True
	name = getTumorName(bam)
	if type(name) != str:
		idx = False
		bam, name = addRG(bam, sid, picard)
		if not bam or not name:
			return None, None
		else:
			idx = True
	if idx == True:
		samIndex(bam)
		return name, bam
