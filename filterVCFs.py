'''This script will filter mutect2 output files.'''

import os
from argparse import ArgumentParser
from datetime import datetime
from glob import glob
from subprocess import Popen
from shlex import split
from commonUtil import *

#-------------------------------Comparison------------------------------------

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
	if "PASS" in outpath:
		typ = "pass"
	else:
		typ = "all"
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
	sa = vcfs[0][vcfs[0].rfind("/")+1:vcfs[0].rfind(".")]
	sb = vcfs[1][vcfs[1].rfind("/")+1:vcfs[1].rfind(".")]
	with open(outpath.replace("_PASS", "") + ".csv", "a") as output:
		output.write(("{},{},{},{},{},{},{}\n").format(typ,sa, sb, a, b, c, sim))
	return 1

def comparePair(outpath, vcfs):
	# Calls gatk and pyvcf to filter and compare given pair of
	done = 0
	# Call bftools on all results
	cmd = ("bcftools isec {} {}").format(vcfs[0], vcfs[1])
	cmd += " -p {}"
	done += intersect(outpath, cmd, vcfs)
	# Call bftools on passes
	outpath += "_PASS"
	done += intersect(outpath, cmd + " -f .,PASS", vcfs)
	return done

def compareVCFs(conf, samples):
	# Compares unfilted vs. filtered results for each combination of pair of samples
	done = 0
	outpath = conf["outpath"] + conf["sample"] + "_"
	with open(outpath + "csv", "w") as output:
		# Initialize summary file and write header
		output.write("Type,SampleA,SampleB,#PrivateA,#PrivateB,#Common,Similarity\n")
	s1, s2 = list(samples.keys())
	# Append filtered smaple name when submitting
	done += comparePair(outpath + sample[s1].ID, samples[s1].Output, samples[s2].Unfiltered)
	done += comparePair(outpath + sample[s2].ID, samples[s2].Output, samples[s1].Unfiltered)
	if done == 4:
		return True
	else:
		return False

#-------------------------------Contamination/Filtering-----------------------

def filterCalls(conf, vcf):
	# Calls gatk to filter mutect calls
	outfile = vcf[:vcf.find(".")] + ".filtered.vcf"
	log = outfile.replace("vcf", "stdout")
	# Assemble command
	if "gatk" in conf.keys():
		cmd = ("java -jar {} FilterMutectCalls ").format(conf["gatk"])
	else:
		cmd = "gatk FilterMutectCalls "
	if "contaminant" in conf.keys():
		cmd += ("-contamination-table {} ").format(conf["contaminant"])
	if "fmo" in conf.keys():
		cmd += " " + conf["fmo"]
	cmd += ("-V {} -O {}").format(vcf, outfile)
	with open(log, "w") as l:
		try:
			fmc = Popen(split(cmd), stdout = l, stderr = l)
			fmc.communicate()
		except:
			print(("\t[Error] Could not call FilterMutectCalls on {}").format(vcf))
			return None
	return bgzip(outfile)

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
	pileup = s.Output.replace(".vcf", ".pileup.table")
	pu += ("-I {} -V {} -O {}").format(s.Input, conf["contaminant"], pileup)
	pu = getOpt(conf, pu)
	plog = pileup.replace("table", "stdout")
	with open(plog, "w") as l:
		try:
			spu = Popen(split(pu), stdout = l, stderr = l)
			spu.communicate()
		except:
			print(("\t[Error] Could not call GetPileupSummaries on {}").format(s.Output))
			return False
	if getStatus(plog) == False:
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
	return getStatus(clog)

#-----------------------------------------------------------------------------

def printError(msg):
	# Prints formatted error message
	print(("\n\t[Error] {}. Skipping.\n").format(msg))

def checkSamples(name, samples):
	# Makes sure two samples have passed mutect
	proceed = True
	msg = ""
	if len(samples.keys()) != 2:
		printError(("Could not find two sample vcfs for {}").format(name))
		proceed =  False
	else:
		# Ensure both have at least passed mutect
		for s in samples.keys():
			if samples[s].Step == "mutect" and samples[s].Status != "complete":
				printError(("{} from {} has not successfully completed mutect.").format(samples[s].ID, name)) 
				proceed = False
			elif samples[s].Step == "filtering" and samples[s].Status != "complete":
				# Re-attempt failed filtering steps
				samples[s].Step = "contamination-estimate"
 	return proceed

def filterSamples():
	# Filters and estimates contamination
	paths = glob(conf["outpath"])
	for p in paths:
		# Iterate through each subdirectory
		log, samples = checkOutput(p)
		conf["log"] = log
		# Get sample name from path
		if sample[:-1] == "/":
			sample = sample[:-1]
		sample = p[p.rfind("/")+1:]
		proceed = checkSamples(sample, samples)
		if proceed == True:
			# Compare output
			print(("\n\tFiltering and comparing VCFs from {}...").format("sample"))
			for s in samples.keys():
				if samples[s].Step == "mutect" and samples[s].Status == "complete":
					samples[s].Step = "contamination-estimate"
					samples[s].Status = "starting"
					if "contaminant" in conf.keys():
						# Estimate contamination
						status = estContamination(conf, samples[s])
						if status == True:
							samples[s].Status = "complete"
						else:
							samples[s].Status = "failed"
					else:
						s.Status = "none"
					appendLog(conf, samples[s])
				if samples[s].Step == "contamination-estimate":
					# Filter vcf if contamination est has been attempted or previous filtering failed
					samples[s].Step = "filtering"
					samples[s].Status = "starting"
					filtered = filterCalls(conf, samples[s].Output)
					if filtered:
						# Record filtered reads
						samples[s].Output = filtered
						samples[s].Status = "complete"
					else:
						samples[s].Status = "failed"
					appendLog(conf, s)
					samples[s].Step = "comparison"
					samples[s].Status = "starting"
			# Comparison
			status = compareVCFs(conf, samples)
			if status == True:
				print(("\tAll files for {} filtered successfully.").format(sample))
				for s in samples.keys():
					samples[s].Status = "complete"
					appendLog(conf, samples[s])
			else:
				print(("\t[Error] Some files from {} failed comparison.").format(sample))
				for s in samples.keys():
					samples[s].Status = "failed"
					appendLog(conf, samples[s])
			
def main():
	starttime = datetime.now()
	parser = ArgumentParser("This script will filter mutect2 output files.")
	parser.add_argument("-c", help = "Path to config file containing reference genome, java jars \
(if using), and mutect options (required; input files are read from sub-directories in output_directory \
and output will be written to same sub-directory).")
	args = parser.parse_args()
	# Load config file and discard batch template
	conf, _ = getConf(args.c)
	done = filterSamples(conf)
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
