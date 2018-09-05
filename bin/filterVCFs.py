'''This script will filter mutect2 output files.'''

import os
from shutil import rmtree
from argparse import ArgumentParser
from datetime import datetime
from glob import glob
from multiprocessing import Pool, cpu_count
from commonUtil import *
from vcfCompare import *
from unixpath import checkDir

def cleanUp(outpath):
	# Removes intermediate files
	paths = glob(outpath + "*")
	for i in paths:
		if os.path.isdir(i):
			if "_unfiltered" in i:
				# Remove intermediate isec results
				rmtree(i)
		else:
			if ".stdout" in i:	
				# Remove logs
				os.remove(i)
			elif ".unfiltered" in i or ".noGermline" in i:
				# Remove unfiltered results
				os.remove(i)

def covB(S):
	# Calls bcf isec and uses output to generate bed files for each comparison
	run = False
	if S.A.Step == "filtering_germline" and S.A.Status == "complete":
		run = True
	elif S.A.Step == "filtering_covB" and S.A.Status != "complete":
		run = True
	if run == True:
		# Update statuses and get output file names and log
		S.updateStatuses("starting", "filtering_covB")
		S.A.Output = samples["A"].Unfiltered.replace(".noGermline", ".covB")
		S.B.Output = samples["B"].Unfiltered.replace(".noGermline", ".covB")
		# Call covB.sh: vcf1 vcf2 outputvcf2 outputvcf1 bam1 bam2 genome gatkjar
		cmd = ("bash covB.sh {} {} {} {} ").format(S.A.Private, samples["B"].Private, samples["B"].Output, samples["A"].Output)
		cmd += ("{} {} {} {}").format(samples["A"].Bam, samples["B"].Bam, conf["ref"], conf["gatk"])
		print(cmd)
		quit()
		res = runProc(cmd, log)
		if res == True:
			# Output names have already been updated
			samples["A"].updateStatus("complete")
			samples["B"].updateStatus("complete")
		else:
			samples["A"].updateStatus("failed")
			samples["B"].updateStatus("failed")
		appendLog(conf, samples["A"])
		appendLog(conf, samples["B"])
	return samples

def rmGermline(conf, sample, outpath):
	# Calls filterMutectCalls and bcftools to remove germline risks
	run = False
	if sample.Step == "mutect" and sample.Status == "complete":
		run = True
	elif sample.Step == "filtering_germline" and sample.Status == "failed":
		run == True
	if run == True:
		sample.updateStatus("starting", "filtering_germline")
		unfiltered, res = filterCalls(conf, sample.Output, "a", outpath)
		if res == True and ".noGermline." in unfiltered:
			# Record unfiltered reads
			sample.updateStatus("complete", outfile = unfiltered, unfilt = True)
		else:
			sample.updateStatus("failed")
		appendLog(conf, sample)	
	return sample	

def filterPair(conf, S):
	# Filters and compares pair of samples
	covb = False
	nan = False
	S.A = rmGermline(S.Conf, S.A, S.Outdir)
	S.B = rmGermline(S.Conf, S.B, S.Outdir)
	print(S.A, S.B)
	quit()
	if S.B.Status == "complete" and S.A.Status == "complete":
		# Add summary to unfiltered log and use output of bcfIsec
		S.compareVCFs()
		# Update statuses
			if status == True:
				S.A.Private = 
				S.updateStatuses("complete", "filtering_isec1", True)
				covb = True
			else:
				S.updateStatuses("failed", "filtering_isec1", True)
	if covb == True:
		S = covB(S)

	if nan == True and conf["cleanup"] == True:
		# Remove intermediary files if indicated and program exited successfully
		cleanUp(variants["outpath"])
	return [nan, variants["ID"]]

#--------------------------------------------I/O------------------------------

def getOutdir(conf, outdir, done, flog, ulog):
	# Reads in dictionary of input samples
	variants = []
	print("\tReading input vcfs...")
	paths = glob(conf["outpath"] + "*")
	for p in paths:
		# Iterate through each subdirectory
		S = Samples()
		S.setLogs(flog, ulog, conf)
		res = S.setSamples(p, outdir, done)
		if res == True:
			variants.append(S)
	return variants

def getComplete(outdir):
	# Makes log file or returns list of completed samples
	summary = self.Outdir + "summary_Filtered.csv"
	ulog = outdir + "summary_Unfiltered.csv"
	for log in [ulog, summary]:
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
	return done, summary, ulog
			
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
		done, flog, ulog = getComplete(args.o)
	else:
		done, flog, ulog = getComplete(conf["outpath"])
	variants = getOutdir(conf, args.o, done, flog, ulog)
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
