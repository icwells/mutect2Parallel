'''This script will filter mutect2 output files.'''

import os
from shutil import copy
from argparse import ArgumentParser
from datetime import datetime
from glob import glob
from functools import partial
from multiprocessing import Pool, cpu_count
from subprocess import Popen
from shlex import split
from commonUtil import *
from bamUtils import *
from unixpath import checkDir

def printError(msg):
	# Prints formatted error message
	print(("\n\t[Error] {}. Skipping.\n").format(msg))

#-------------------------------Comparison------------------------------------

def comparePair(outpath, vcfs):
	# Calls gatk and pyvcf to filter and compare given pair of
	for i in range(len(vcfs)):
		if vcfs[i] and os.path.isfile(vcfs[i]):
			# Make sure files are bgzipped
			vcfs[i] = tabix(vcfs[i])
		elif ".gz" not in vcfs[i] and os.path.isfile(vcfs[i] + ".gz"):
			vcfs[i] += ".gz"
		else:
			printError(("Cannot find {}").format(vcfs[i]))
			return None
	# Call bftools on all results
	ret = bcfIsec(outpath, vcfs)
	return ret

def compareVCFs(conf, log, name, samples):
	# Compares unfilted vs. passed results for each combination of pair of samples
	cont = False
	ret = False
	outpath = conf["outpath"] + name + "/"
	print("\tComparing samples...")
	s1, s2 = list(samples.keys())
	# Append filtered sample name when submitting
	aout = outpath + "_" + samples[s1].ID
	a = comparePair(aout, [samples[s1].Output, samples[s2].Unfiltered])
	if a:
		bout = outpath + "_" + samples[s2].ID
		b = comparePair(bout, [samples[s2].Output, samples[s1].Unfiltered])
		if b:	
			# Only continue if both pass
			cont = True
	if cont == True:
		# Merge common variants and get total and similarity
		common = bcfMerge(outpath, [acom, bcom])
		if common:
			c = getTotal(common)
			try:
				sim = c/(a+b+c)
			except ZeroDivisionError:
				sim = 0.0
			with open(log, "a") as out:
				out.write(("{},{},{},{},{},{},{:.2%}\n").format(name, samples[s1].ID, samples[s2].ID, a, b, c, sim))
			ret = True
	return ret

#-------------------------------Filtering-------------------------------------

def bcftoolsFilter(vcf):
	# Calls bcftools to filter calls before calling isec
	outfile = vcf.replace("unfiltered.vcf.gz", "passed.vcf")
	cmd = ('bcftools filter -e "FILTER=\'germline_risk\'" -o {} {}').format(outfile, vcf)
	try:
		fmc = Popen(split(cmd))
		fmc.communicate()
		return tabix(outfile)
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
		return tabix(outfile)
	else:
		return None

#-----------------------------------------------------------------------------

def filterPair(conf, log, variants):
	# Filters and compares pair of samples
	compare = False
	conf["log"] = variants["log"]
	print(("\n\tFiltering and comparing VCFs from {}...").format(sample))
	samples = variants["samples"]
	for s in samples.keys():
		if samples[s].Step == "mutect" and samples[s].Status == "complete" and not conf["summary"]:
			print(("\tFiltering {}...").format(samples[s].ID))
			# FilterMutectCalls
			samples[s].Step = "filtering"
			samples[s].Status = "starting"
			unfiltered = filterCalls(conf, samples[s].Output, variants["outpath"])
			if unfiltered:
				# Record unfiltered reads
				samples[s].Output = unfiltered
				samples[s].Unfiltered = unfiltered
			else:
				samples[s].Status = "failed"
				appendLog(conf, samples[s])
			if samples[s].Status != "failed":
				# bcftools filter
				passed = bcftoolsFilter(samples[s].Output)
				if passed:
					samples[s].Output = passed
					samples[s].Status = "complete"
					appendLog(conf, samples[s])
					samples[s].Step = "comparison"
					samples[s].Status = "starting"
					compare = True
				else:
					samples[s].Status = "failed"
					appendLog(conf, samples[s])
	if compare == True:
		# Comparison
		status = compareVCFs(conf, log, variants["ID"], samples)
		if status == True:
			for s in samples.keys():
				samples[s].Status = "complete"
				appendLog(conf, samples[s])
		else:
			compare = False
	if compare == False:
		# Call if filtering/comparison failed
		for s in samples.keys():
			samples[s].Status = "failed"
			appendLog(conf, samples[s])
	return [compare, variants["ID"]]

#--------------------------------------------I/O------------------------------

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

def getOutdir(conf, outdir, done):
	# Reads in dictionary of input samples
	variants = []
	print("\tReading input vcfs...")
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
			# Get sample name from path sans trailing slash
			sid = p[:-1]
			sid = sid[sid.rfind("/")+1:]
			if sid not in done:
				# Proceed if sample not done and log is present
				log, s = checkOutput(p, False)
				proceed = checkSamples(sid, s)
				if proceed == True:
					sampleout = checkDir(conf["outpath"] + sid, True)
					if not os.path.isfile(sampleout + "mutectLog.txt"):
						copy(log, sampleout + "mutectLog.txt")
					# {ID, samples, log file, outdir}
					d = {"ID": sid, "samples": s, "log": sampleout + "mutectLog.txt", "outpath": sampleout}
					variants.append(d)
	return variants

def getComplete(log):
	# Makes log file or returns list of completed samples
	first = True
	done = []
	if not os.path.isfile(log):
		with open(log, "w") as out:
			# Initialize summary file and write header
			out.write("ID,SampleA,SampleB,#PrivateA,#PrivateB,#Common,%Similarity\n")
	else:
		with open(log, "r") as f:
			for line in f:
				if first == False:
					done.append(line.split(",")[0])
				else:
					first = False
	return done
			
def main():
	starttime = datetime.now()
	parser = ArgumentParser("This script will filter mutect2 output files.")
	parser.add_argument("-t", type = int, default = 1, 
help = "Number of threads.")
	parser.add_argument("-c", help = "Path to config file containing reference genome, java jars \
(if using), and mutect options (required; input files are read from sub-directories in output_directory \
and output will be written to same sub-directory).")
	parser.add_argument("-o", help = "Output directory (if different from directory in config file).")
	args = parser.parse_args()
	if args.t > cpu_count():
		args.t = cpu_count()
	# Load config file and discard batch template
	conf, _ = getConf(args.c)
	log = conf["outpath"] + "summary.csv"
	done = getComplete(log)
	variants = getOutdir(conf, args.o, done)
	pool = Pool(processes = args.t)
	func = partial(filterPair, conf, log)
	print(("\tComparing samples with {} threads...").format(args.t))
	for x in pool.imap_unordered(func, variants):
		if x[0] == False:
			print(("\t[Error] Some files from {} failed comparison.").format(x[1]))
		else:		
			print(("\tAll comparisons for {} run successfully.").format(x[1]))
	pool.close()
	pool.join()
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
