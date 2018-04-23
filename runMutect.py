'''This script contains functions for calling mutect2 and filtering its output'''

import os
from subprocess import Popen
from shlex import split
from bamUtil import *
from runPair import Sample

#-------------------------------Filtering-------------------------------------

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
	done += intersect(outpath, cmd + " -f .,PASS", vcfs)
	if done == 2:
		ret = True
	return ret

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
			return None
	return bgzip(outfile)

#-------------------------------Contamination---------------------------------

def estContamination(conf, s):
	# Calls gatk to get pileup summary and estimate contamination
	print(("\tEstimating contamination in {}...").format(s.ID))
	if "gatk" in conf.keys():
		# Format command for calling gatk jar
		pu = ("java -jar {} GetPileupSummaries -R {} ").format(conf["gatk"], conf["reference"])
		cc = ("java -jar {} CalculateContamination ").format(conf["gatk"])
	else:
		# Format command for colling gatk from path
		pu = ("gatk GetPileupSummaries -R {} ").format(conf["reference"])
		cc = "gatk CalculateContamination "
	# Get pileup summary
	pileup = s.Output.replace(".vcf", "pileup.table")
	pu += ("-I {} -V {} -O {}").format(s.Output, conf["contaminant"], pileup)
	pu = getOpt(conf, pu)
	plog = pileup.replace("table", "stdout")
	with open(plog, "w") as l:
		try:
			spu = Popen(split(pu), stdout = l, stderr = l)
			spu.communicate()
		except:
			print(("\t[Error] Could not call GetPileupSummaries on {}").format(s.Output))
			return False
	# Get contamination estimate
	cest = pileup.replace("pileup", "contamination")
	clog = cest.replace("table", "stdout")
	cc += ("-I {} -O {}").format(pileup, cest)
	with open(clog, "w") as l:
		try:
			ccont = Popen(split(cc), stdout = l, stderr = l)
			ccont.communicate()
		except:
			print(("\t[Error] Could not call CalculateContamination on {}").format(s.Output))
			return False
	return True

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
		cmd += (" --normal_panel {}").format(conf["pon"])
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
	if s.Status == "starting":
		s = submitSample(infile, conf, s, name)
	if s.Status == "mutect":
		if "contaminant" in conf.keys():
			status = estContamination(conf, s)
			if status == True:
				s.Status = "contamination-estimate:complete"
			else:
				s.Status = "failed:estimating-contamination"
		else:
			s.Status = "contamination-estimate:none"
		appendLog(conf, s)
	if "contamination-estimate" in s.Status:
		# Assemble command
		if "gatk" in conf.keys():
			cmd = ("java -jar {} FilterMutectCalls ").format(conf["gatk"])
		else:
			cmd = "gatk FilterMutectCalls "
		# Filter vcf
		filtered = filterCalls(cmd, s.Output)
		if filtered:
			# Record filtered reads
			s.Output = filtered
			s.Status = "filtered"
		else:
			s.Status = "failed-filtering"
		appendLog(conf, s)
	return s
