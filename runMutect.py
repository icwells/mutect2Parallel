'''This script contains functions for calling mutect2 and filtering its output'''

import os
from subprocess import Popen
from shlex import split
from bamUtil import *
from runPair import Sample

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
