'''This script will filter mutect2 output files.'''

import os
from shutil import copy
from argparse import ArgumentParser
from datetime import datetime
from glob import glob
from subprocess import Popen
from shlex import split
from commonUtil import *
from unixpath import checkDir

def printError(msg):
	# Prints formatted error message
	print(("\n\t[Error] {}. Skipping.\n").format(msg))

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

def intersect(summary, log, outpath, cmd, vcfs):
	# Calls bcftools to get intersecting rows and summarizes output
	summarize = False
	if summary == False:
		try:
			bcf = Popen(split(cmd))
			bcf.communicate()
			summarize = True
		except:
			print(("\t[Error] Could not call bcftools isec with {}").format(cmd))
			return 0
	else:
		summarize = True
	if summarize == True:
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
		with open(log, "a") as output:
			output.write(("{},{},{},{},{},{}\n").format(sa, sb, a, b, c, sim))
		return 1

def comparePair(summary, log, outpath, vcfs):
	# Calls gatk and pyvcf to filter and compare given pair of
	for i in range(len(vcfs)):
		if vcfs[i] and os.path.isfile(vcfs[i]):
			# Make sure files are bgzipped
			vcfs[i] = bgzip(vcfs[i])
		elif ".gz" not in vcfs[i] and os.path.isfile(vcfs[i] + ".gz"):
			vcfs[i] += ".gz"
		else:
			printError(("Cannot find {}").format(vcfs[i]))
			return False
	# Call bftools on all results
	cmd = ("bcftools isec {} {} -p {}").format(vcfs[0], vcfs[1], outpath)
	done = intersect(summary, log, outpath, cmd, vcfs)
	return done

def compareVCFs(conf, name, samples):
	# Compares unfilted vs. passed results for each combination of pair of samples
	done = 0
	outpath = conf["outpath"] + name + "/" + name
	log = outpath + ".csv"
	print("\tComparing samples...")
	with open(log, "w") as output:
		# Initialize summary file and write header
		output.write("SampleA,SampleB,#PrivateA,#PrivateB,#Common,%Similarity\n")
	s1, s2 = list(samples.keys())
	# Append filtered smaple name when submitting
	done += comparePair(conf["summary"], log, outpath + "_" + samples[s1].ID, [samples[s1].Output, samples[s2].Unfiltered])
	done += comparePair(conf["summary"], log, outpath + "_" + samples[s2].ID, [samples[s2].Output, samples[s1].Unfiltered])
	if done == 2:
		return True
	else:
		return False

#-------------------------------Filtering-------------------------------------

def bcftoolsFilter(vcf):
	# Calls bcftools to filter calls before calling isec
	outfile = vcf.replace("unfiltered.vcf.gz", "passed.vcf")
	opt = '"."|","|"PASS"'
	cmd = ("bcftools filter -i 'FILTER={}' -o {} {}").format(opt, outfile, vcf)
	try:
		fmc = Popen(split(cmd))
		fmc.communicate()
		return bgzip(outfile)
	except:
		print(("\t[Error] Could not call bcftools filter on {}").format(vcf))
		return None

def filterCalls(conf, vcf, outdir = None):
	# Calls gatk to filter mutect calls to remove germline variants
	if outdir:
		outfile = outdir + vcf[vcf.rfind("/")+1:vcf.find(".")] + ".unfiltered.vcf"
	else:
		outfile = vcf[:vcf.find(".")] + ".unfiltered.vcf"
	log = outfile.replace("vcf", "stdout")
	# Assemble command
	if "gatk" in conf.keys():
		cmd = ("java -jar {} FilterMutectCalls ").format(conf["gatk"])
	else:
		cmd = "gatk FilterMutectCalls "
	'''if "contaminant" in conf.keys():
		cmd += ("-contamination-table {} ").format(conf["contaminant"])
	if "fmo" in conf.keys():
		cmd += " " + conf["fmo"]'''
	cmd += ("-V {} -O {}").format(vcf, outfile)
	with open(log, "w") as l:
		try:
			fmc = Popen(split(cmd), stdout = l, stderr = l)
			fmc.communicate()
		except:
			print(("\t[Error] Could not call FilterMutectCalls on {}").format(vcf))
	if getStatus(log) == True:
		return bgzip(outfile)
	else:
		return None

#-----------------------------------------------------------------------------

def filterSamples(conf, samples):
	# Filters and compares all vcfs in each subdirectory
	for sample in samples.keys():
		# Compare output
		pair = []
		compare = False
		conf["log"] = samples[sample]["log"]
		print(("\n\tFiltering and comparing VCFs from {}...").format(sample))
		for i in samples[sample].keys():
			if i != "log" and i != "outpath":
				pair.append(samples[sample][i])
		for s in pair:
			if s.Step == "mutect" and s.Status == "complete" and not conf["summary"]:
				print(("\tFiltering {}...").format(s.ID))
				# FilterMutectCalls
				s.Step = "filtering"
				s.Status = "starting"
				unfiltered = filterCalls(conf, s.Output, samples[sample]["outpath"])
				if unfiltered:
					# Record unfiltered reads
					s.Output = unfiltered
					s.Unfiltered = unfiltered
				else:
					s.Status = "failed"
					appendLog(conf, s)
				if s.Status != "failed":
					# bcftools filter
					passed = bcftoolsFilter(s.Output)
					if passed:
						s.Output = passed
						s.Status = "complete"
						appendLog(conf, s)
						s.Step = "comparison"
						s.Status = "starting"
						compare = True
					else:
						s.Status = "failed"
						appendLog(conf, s)
		if compare == True:
			# Comparison
			status = compareVCFs(conf, sample, samples[sample])
			if status == True:
				print(("\tAll comparisons for {} run successfully.").format(sample))
				for s in pair:
					s.Status = "complete"
					appendLog(conf, s)
			else:
				compare = False
		if compare == False:
			# Call if filtering/comparison failed
			print(("\t[Error] Some files from {} failed comparison.").format(sample))
			for s in pair:
				s.Status = "failed"
				appendLog(conf, s)

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
			if not os.path.isfile(samples[s].Output):
				if os.path.isfile(samples[s].Output + ".gz"):
					samples[s].Output = samples[s].Output + ".gz"
				else:
					printError(("Cannot find {} from {}").format(samples[s].Output, name)) 
					proceed = False	
			if samples[s].Step == "mutect":
				if samples[s].Status != "complete":
					printError(("{} from {} has not successfully completed mutect.").format(samples[s].ID, name)) 
					proceed = False
				elif ".vcf" not in samples[s].Unfiltered:
					printError(("Cannot find mutect output for {} from {}.").format(samples[s].ID, name)) 
					proceed = False					
			elif samples[s].Step == "filtering" and samples[s].Status != "complete":
				# Re-attempt failed filtering steps
				samples[s].Step = "mutect"
	return proceed

def getOutdir(conf, outdir):
	# Reads in dictionary of input samples
	samples = {}
	print("\n\tReading input vcfs...")
	paths = glob(conf["outpath"] + "*")
	if outdir and outdir != conf["outpath"]:
		outdir = checkDir(outdir)
		# Reassign ouutpath after getting input directories
		conf["outpath"] = outdir
	for p in paths:
		# Iterate through each subdirectory
		if p[:-1] != "/":
			p += "/"
		if os.path.isfile(p + "mutectLog.txt"):
			# Proceed if log is present
			log, s = checkOutput(p, False)
			# Get sample name from path sans trailing slash
			sid = p[:-1]
			sid = sid[sid.rfind("/")+1:]
			proceed = checkSamples(sid, s)
			if proceed == True:
				sampleout = checkDir(conf["outpath"] + sid, True)
				if not os.path.isfile(sampleout + "mutectLog.txt"):
					copy(log, sampleout + "mutectLog.txt")
				# {ID: {sample 1, sample 2, log file, outdir}
				samples[sid] = s
				samples[sid]["log"] = sampleout + "mutectLog.txt"
				samples[sid]["outpath"] = sampleout
	return samples
			
def main():
	starttime = datetime.now()
	parser = ArgumentParser("This script will filter mutect2 output files.")
	parser.add_argument("-c", help = "Path to config file containing reference genome, java jars \
(if using), and mutect options (required; input files are read from sub-directories in output_directory \
and output will be written to same sub-directory).")
	parser.add_argument("-o", help = "Output directory (if different from directory in config file).")
	parser.add_argument("--summarize", action = "store_true", default = False,
help = "Skips to summarize step (all other steps must be completed. overwrites existing summary files).")
	args = parser.parse_args()
	# Load config file and discard batch template
	conf, _ = getConf(args.c)
	conf["summary"] = args.summarize
	samples = getOutdir(conf, args.o)
	done = filterSamples(conf, samples)
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
