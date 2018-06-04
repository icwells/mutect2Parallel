'''This script contains functions for adding readgroups to bam files, as well as indexing and extracting readgroups'''

import os
import pysam
from subprocess import Popen
from shlex import split

class Sample():
	# Stores data for managing sample progress
	def __init__(self, name, step, status, outfile):
		self.ID = name
		self.Step = step
		self.Status = "starting"
		if "failed" not in status:
			self.Status = status
		self.Output = outfile
		self.Bam = None
		self.Input = None

	def __newStatus__(self, status):
		# Updates self.Status
		if "failed" not in status:
			self.Status = status
		else:
			self.Status = "starting"

	def __update__(self, name, step, status, outfile):
		# Sorts and updates entry with additional status update
		if step == "filtering" and self.Step == "mutect":
			self.Step = step
			self.Output = outfile
			self.__newStatus__(status)
		elif step == "mutect" and self.Step == "starting":
			self.Status = mostRecent
			self.Output = outfile
			self.__newStatus__(status)

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
			l.write(("{}\t{}\t{}\n").format(s.ID, s.Status, s.Output))

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
					if len(line) == 4:
						if line[0] in done.keys():
							done[line[0]].__update__(line[0], line[1], line[2], line[3])
						else:
							done[line[0]] = Sample(line[0], line[1], line[2], line[3])
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

#-------------------------------bamfunctions----------------------------------------

def bgzip(vcf):
	# bgzip compresses filtered vcf files
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	H_SL-20-LK38-020-node-A1	H_SL-20-LK40-020-B3
	if os.path.isfile(vcf + ".gz"):
		gz = vcf + ".gz"
	else:	
		gz = pysam.tabix_index(vcf, seq_col=0, start_col=1, end_col=1)
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
