'''This script contains functions for adding readgroups to bam files, as well as indexing and extracting readgroups'''

import os
import gzip
import pysam
from sys import stderr
from subprocess import Popen
from shlex import split
from unixpath import *
from sample import Sample

def runProc(cmd, log = None):
	# Wraps call to Popen, writes stdout/stdout err to log/devnull, returns True if no errors
	if not log:
		log = os.devnull
	with open(log, "w") as out:
		try:
			call = Popen(split(cmd), stdout = out, stderr = out)
			call.wait()
			if call.returncode is not None:
				return True
		except:
			s = cmd.split()
			proc = s[0]
			if "-" not in s[1]:
				# Get subcommand if present
				proc += " " + s[1]
			elif proc == "java":
				# Replace call to jar with name of jar
				proc = getFileName(s[2])
			print(("\t[Warning] Could not call {}").format(proc), file=stderr)
			return False

def tabix(vcf, force = False, keep = False):
	# tabix index and bgzips vcf files
	if force == False and os.path.isfile(vcf + ".gz"):
		gz = vcf + ".gz"
	else:
		vcf = checkGZ(vcf)
		try:
			gz = pysam.tabix_index(vcf, seq_col=0, start_col=1, end_col=1, force=True, keep_original=keep)
		except OSError:
			print(("\t[Warning] Could not index {}.").format(vcf), file=stderr)
			gz = None
	return gz

def bcfMerge(outfile, com):
	# Calls bftools concat on input files
	cmd = ("bcftools merge --force-samples -O v -o {} {} {}").format(outfile, com[0], com[1])
	res = runProc(cmd)
	if res == False:
		outfile = None
	return outfile	

def bcfSort(infile):
	# Call bcftools sort
	if not infile:
		return None
	infile = checkGZ(infile)
	outfile = infile.replace(".vcf", ".sorted.vcf")
	cmd = ("bcftools sort -o {} {}").format(outfile, infile)
	res = runProc(cmd)
	if res == False or os.path.isfile(outfile) == False:
		return None
	return tabix(outfile, True)

def getTotal(vcf):
	# Returns total number of content lines from vcf
	if os.path.isfile(vcf):
		if getExt(vcf) == "gz":
			try:
				count = 0
				with gzip.open(vcf, "rb") as f:
					tag = ("#").encode()
					for line in f:
						if tag not in line:
							count += 1
			except (OSError, UnicodeDecodeError):
				print(("\t[Error] Could not read {}").format(vcf), file=stderr)
				count = None
		else:
			count = 0
			with open(vcf, "r") as f:
				for line in f:
					if line[0] != "#":
						count += 1
	else:
		print(("\t[Error] Could not find {}").format(vcf), file=stderr)
		count = None
	return count

def bcfIsec(outpath, vcfs):
	# Calls bcftools to get intersecting rows and returns number of private A
	a = None
	for i in range(len(vcfs)):
		# Make sure there is an up-to-date index file
		vcfs[i] = tabix(vcfs[i], force = True)
	if None not in vcfs:
		cmd = ("bcftools isec {} {} -p {}").format(vcfs[0], vcfs[1], outpath)
		res = runProc(cmd)
		if res == True:
			# Number of unique variants to each sample and number of shared
			a = getTotal(outpath + "/0000.vcf")
	return a

#-----------------------------------------------------------------------------

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
	print(("\tAdding read groups to {}").format(bam))
	res = runProc(cmd)
	if res == True and outfile:
		name = getTumorName(outfile)
		return outfile, name
	else:
		return None, None

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
#-------------------------------commonfunctions----------------------------------------

def checkGZ(f):
	# Checks for unzipped/gzipped file
	if not os.path.isfile(f):
		if getExt(f) != "gz":
			if os.path.isfile(f + ".gz"):
				# Change infile if file is properly named
				f += ".gz"
		else:
			if os.path.isfile(f + ".gz"):
				# Fix filename on system
				os.rename(f + ".gz", f)
			elif os.path.isfile(f.replace(".gz", "")):
				f = f.replace(".gz", "")
	elif f.count(".gz") > 1:
			old = f
			f = f.replace(".gz.gz", ".gz")
			os.rename(old, f)
	return f

def printError(msg):
	# Prints formatted error message
	print(("\n\t[Error] {}. Skipping.\n").format(msg), file=stderr)

def getDelim(line):
	# Returns delimiter
	for i in ["\t", ",", " "]:
		if i in line:
			return i
	print("\n\t[Error] Cannot determine delimeter. Check file formatting. Exiting.\n", file=stderr)
	quit()

def getStatus(log):
	# Returns true if tool returns success
	status = False
	ret = False
	with open(log, "r") as l:
		for line in l:
			if "error" in line.lower():
				break
			if "Tool returned:" in line:
				status = True
			elif status == True:
				# Get exit status
				if "SUCCESS" in line:
					ret = True
				break
	return ret

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

def checkOutput(outdir, normal = None, prnt = True):
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
			if normal:
				f.write("N\t{}\tnormal\tcomplete\t{}\n".format(getFileName(normal), normal))
	return log, done

def configEntry(conf, arg, key):
	# Returns dict with updated arg entry
	if not arg:
		print(("\n\t[Error] Please specify {}. Exiting.\n").format(arg), file=stderr)
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
			conf["outpath"] = val
		elif target == "GATK_jar":
			conf["gatk"] = val
			checkFile(conf["gatk"])
		elif target == "Picard_jar":
			conf["picard"] = val
			checkFile(conf["picard"])
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
		elif target == "min_covA":
			conf["min_covA"] = int(val)
		elif target == "min_reads_strand":
			conf["min_reads_strand"] = int(val)
		elif target == "min_covB":
			conf["min_covB"] = int(val)
		elif target == "max_altB":
			conf["max_altB"] = int(val)
		elif target == "max_prop_altB":
			conf["max_prop_altB"] = float(val)
		elif target == "max_covN":
			conf["min_covN"] = int(val)
		elif target == "min_freq_altN":
			conf["max_freq_altN"] = float(val)
		elif target == "min_reads_altN":
			conf["max_reads_altN"] = int(val)
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
	checkFile(conf["ref"])
	conf["outpath"] = checkDir(conf["outpath"], True)
	checkFile(conf["gatk"])
	if "germline" in conf.keys() and "af" not in conf.keys():
		print("\n\t[Error] Please supply an allele frequency when using a germline estimate. Exiting.\n", file=stderr)
		quit()	
	return conf, batch
