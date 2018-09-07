'''This script will filter mutect2 output files.'''

import os
from shutil import rmtree
from argparse import ArgumentParser
from datetime import datetime
from glob import glob
from multiprocessing import Pool, cpu_count
from commonUtil import *
from samples import *
from sample import *
from unixpath import checkDir

def cleanUp(outpath):
	# Removes intermediate files
	paths = glob(outpath + "*")
	for i in paths:
		if os.path.isdir(i):
			if "_nab" not in i:
				# Remove intermediate isec results
				rmtree(i)
		else:
			if ".stdout" in i:	
				# Remove logs
				os.remove(i)
			elif ".unfiltered" in i or ".noGermline" in i or ".covB" in i:
				# Remove unfiltered results
				os.remove(i)
			elif "common" in i and "_nab" not in i:
				# Remove intermediary isec unions
				os.remove(i)
			elif "normalVariants" in i:
				# Remove variants extracted from normal bam
				os.remove(i)

def filterPair(S):
	# Filters and compares pair of samples
	covb = False
	nab = False
	S.rmGermline()
	if S.B.Status == "complete" and S.A.Status == "complete":
		# Add summary to unfiltered log and use output of bcfIsec
		S.compareVCFs("a")
		# Update statuses
		S.updateStatuses("complete", "isec1", True)
		covb = True
	if covb == True:
		# Filter for coverage in paired tumor sample
		S.covB()
		if S.B.Status == "complete" and S.A.Status == "complete":
			S.filterForB()
			if S.B.Status == "complete" and S.A.Status == "complete":
				# Add summary to log b and use output of bcfIsec
				S.compareVCFs("b")
				# Update statuses
				S.updateStatuses("complete", "isec2", True)
				nab = True
	if nab == True: 
		# Filter for coverage in normal file
		S.NAB()
		'''if S.B.Status == "complete" and S.A.Status == "complete":
			s.filterForN()
			if S.B.Status == "complete" and S.A.Status == "complete":
				# Add isec results to summary and update log
				S.compareVCFs("n")
				# Update statuses
				S.updateStatuses("complete", "isec3", True)
				if conf["cleanup"] == True:
					# Remove intermediary files if indicated and program exited successfully
					cleanUp(variants["outpath"])
		else:
			nab = False'''
	return [nab, S.ID]

#--------------------------------------------I/O------------------------------

def getOutdir(conf, outdir, done, flog, blog, ulog):
	# Reads in dictionary of input samples
	variants = []
	print("\tReading input vcfs...")
	paths = glob(conf["outpath"] + "*")
	for p in paths:
		if os.path.isfile(p) == False:
			# Iterate through each subdirectory
			S = Samples()
			S.setLogs(flog, ulog, blog, conf)
			res = S.setSamples(p, outdir, done)
			if res == True:
				variants.append(S)
	return variants

def getComplete(outdir):
	# Makes log file or returns list of completed samples
	summary = outdir + "summary_NAB.csv"
	ulog = outdir + "summary_Unfiltered.csv"
	blog = outdir + "summary_covB.csv"
	for log in [ulog, blog, summary]:
		# Check summary last to keep final output
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
	return done, summary, blog, ulog

def checkBin():
	# Makes sure heterAnalyzer and bash scripts are present in working directory
	for idx,i in enumerate(["heterAnalyzer", "covB.sh", "covN.sh"]):
		if not os.path.isfile(i):
			print(("\n\t[Error] {} not found.").format(i))
			if idx == 0:
				print(("\tRun install.sh to install {}.").format(i))
			print("\tExiting.\n")
			quit()
			
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
	checkBin()
	if args.t > cpu_count():
		args.t = cpu_count()
	# Load config file and discard batch template
	conf, _ = getConf(args.c)
	conf["cleanup"] = args.cleanup
	if args.o:
		args.o = checkDir(args.o, True)
		done, flog, blog, ulog = getComplete(args.o)
	else:
		args.o = conf["outpath"]
		done, flog, blog, ulog = getComplete(conf["outpath"])
	variants = getOutdir(conf, args.o, done, flog, blog, ulog)
	l = len(variants)
	pool = Pool(processes = args.t)
	print(("\tComparing samples from {} sets with {} threads...").format(l, args.t))
	for x in pool.imap_unordered(filterPair, variants):
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
