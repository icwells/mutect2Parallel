'''This script will filter mutect2 output files.'''

import os
from shutil import copy
from argparse import ArgumentParser
from datetime import datetime
from glob import glob
from functools import partial
from multiprocessing import Pool, cpu_count
from commonUtil import *
from vcfCompare import *
from unixpath import checkDir

def rmGermline(conf, sample):
	# Calls filterMutectCalls to remove germline risks
	if sample.Step == "mutect" and sample.Status == "complete":
		sample.updateStatus("starting", "filtering_germline")
		infile = sample.Output
		unfiltered = filterCalls(conf, infile, False, variants["outpath"])
		if unfiltered:
			# Record unfiltered reads
			sample.updateStatus("complete", "filtering_germline", unfiltered, True)
		else:
			sample.updateStatus("failed")
			appendLog(conf, samples[s])	
	return sample	

def filterPair(conf, log, variants):
	# Filters and compares pair of samples
	compare = False
	conf["log"] = variants["log"]
	samples = variants["samples"]
	for s in samples.keys():
		samples[s] = rmGermline(conf, samples[s])
	


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
		# Reassign outpath after getting input directories
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
	return variants, conf

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
	parser.add_argument("--cleanup", action = "store_true", default = False,
help = "Remove intermediary files (default is to keep them).")
	args = parser.parse_args()
	if args.t > cpu_count():
		args.t = cpu_count()
	# Load config file and discard batch template
	conf, _ = getConf(args.c)
	conf["cleanup"] = args.cleanup
	if args.o:
		args.o = checkDir(args.o, True)
		log = args.o + "summary.csv"
	else:
		log = conf["outpath"] + "summary.csv"
	done = getComplete(log)
	variants, conf = getOutdir(conf, args.o, done)
	l = len(variants)
	pool = Pool(processes = args.t)
	func = partial(filterPair, conf, log)
	print(("\tComparing samples from {} sets with {} threads...").format(l, args.t))
	for x in pool.imap_unordered(func, variants):
		l -= 1
		if x[0] == False:
			print(("\t[Error] Some files from {} failed comparison.").format(x[1]))
		else:		
			print(("\tAll comparisons for {} run successfully. {} samples remaining.").format(x[1], l))
	pool.close()
	pool.join()
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
