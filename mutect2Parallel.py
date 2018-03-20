'''This script will call MuTect2 on a given list of input files'''

import argparse
import os
from datetime import datetime
from functools import partial
from multiprocessing import Pool, cpu_count
from subprocess import Popen
from shlex import split
from bamUtil import *

def callMutect(cmd, path, n, t):
	# Calls Mutect with given root command and files
	nid = n[n.rfind("/")+1:].replace(".bam", "")
	tid = t[t.rfind("/")+1:].replace(".bam", "")
	outfile = ("{}-{}.vcf").format(path, tid)
	tumorname, t = checkRG(t, tid)
	cmd += ("--tumor-sample {} -I {} -I {} --output {}").format(tumorname, t, n, outfile)
	print(("\tComparing {} and {}...").format(nid, tid))
	with open(outfile.replace("vcf", "stdout"), "w") as dn:
		try:
			mt = Popen(split(cmd), stdout = dn, stderr = dn)	
			mt.communicate()
			print(("\tFinished comparing {} and {}...").format(nid, tid))
			return outfile
		except:
			print(("\t[Error] Could not call MuTect on {} and {}.").format(nid, tid))
			return ""

def submitFiles(conf, outfiles, outdir, sample):
	# Calls MuTect2 serially over input files
	if sample[4] < 3:
		# Proceed if at least one combination has not been run
		if sample[0] in outfiles.keys():
			vcfs = outfiles[sample[0]]
		else:
			vcfs = []
		print(("\n\tProcessing {}...").format(sample[0]))
		# Assemble command
		if conf["jar"] == False:
			# Format command for colling gatk from path
			cmd = ("gatk Mutect2 -R {} ").format(conf["ref"])
			_, control = checkRG(sample[1], sample[0])
		else:
			# Format command for calling gatk jar
			cmd = ("java -jar {} Mutect2 -R {} ").format(conf["gatk"], conf["ref"])
			_, control = checkRG(sample[1], sample[0], conf["picard"])
		if sample[4] == 1:
			s = [sample[3]]
		elif sample[4] == 2:
			s = [sample[2]]
		else:
			s = [sample[2], sample[3]]
		# Call for each remaining combination of control-tumor
		for j in s:
			res = callMutect(cmd, outdir + sample[0], control, j)
			if res:
				# Record finished samples
				vcfs.append(res)
				with open(conf["log"], "a") as l:
					l.write(("{}\tMutect\n").format(res[res.rfind("/")+1:res.rfind(".")]))
	'''elif sample[4] == 3:
		# Compare output
		status = compareVCFs(vcfs)
		if status = True:
			# Record finished samples
			with open(conf["log"], "a") as l:
				l.write(("{}\tcomparison\n").format(sample[0]))'''
	return sample[0]

def getManifest(done, infile):
	# Returns dict of input files
	files = []
	first = True
	print("\tReading input file...")
	with open(infile, "r") as f:
		for line in f:
			a = 0
			b = 0
			if first == True:
				# Determine delimiter
				for i in [" ", "\t", ","]:
					if i in line:
						delim = i
						break
				first = False
			s = line.strip().split(delim)
			if s[0] in done.keys():
				if done[s[0]][-1] == "comparison":
					# Skip completed files
					pass
				else:
					# [ID, Normal, A, B, status]
					if s[2] in done[s[0]]:
						a = 1
					if s[3] in done[s[0]]:
						b = 2
			# Statuses: 0=none, 1=a_done, 2=b_done, 3=both
			s.append(a+b)
			files.append(s)
	for i in files:
		for j in i[1:4]:
			if not os.path.isfile(j):
				print(("\n\t[Error] Input file {} not found. Exiting.\n").format(j))
				quit()
	print(("\tFound entries for {} samples.").format(len(files)))
	return files

def checkOutput(outdir):
	# Checks for output log file and reads if present
	first = True
	done = {}
	outfiles = {}
	log = outdir + "mutectLog.txt"
	print("\tChecking for previous output...")
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	if os.path.isfile(log):
		with open(log, "r") as f:
			for line in f:
				if first == False:
					line = line.strip().split("\t")
					splt = line[0].split("-")
					if splt[0] in done.keys():
						# Record vcf location
						outfiles[splt[0]].append(splt[2])
						if line[1] == "comparison":
							# Record last stage completed (comparison/mutect)
							done[splt[0]] = [done[splt[0]][0], splt[1], line[1]]
						else:
							done[splt[0]] = [done[splt[0]][0], splt[1], done[splt[0]][1]]
					else:
						# {sample#: [filename1, filename1, mostRecent]}
						done[splt[0]] = [splt[1], line[1]]
						# {sample#: [vcf1, vcf2]}
						outfiles[splt[0]] = [splt[2]]
				else:
					# Skip header
					first = False
		print(("\tFound {} completed samples.").format(len(done)))
	else:
		with open(log, "w") as f:
			# Initialize log file
			f.write("Sample\tStatus\tOutput\n")
	return log, done, outfiles

def checkReferences(conf):
	# Ensures fasta and vcf index and dict files are present
	ref = conf["ref"]
	if not os.path.isfile(ref + ".fai"):
		getFastaIndex(ref)
	fdict = ref.replace(".fa", ".dict")
	with open(os.devnull, "w") as dn:
		if not os.path.isfile(fdict):
			# Call CreateSequenceDictionary
			print("\tGenerating fasta dictionary...\n")
			if conf["jar"] == True:
				cmd = ("java -jar {} ").format(conf["picard"])
			else:
				cmd = "picard "
			fd = Popen(split(("{} CreateSequenceDictionary R= {} O= {}").format(cmd, ref, fdict)), stdout=dn, stderr=dn)
			fd.communicate()

def getConf(cpu, ref, infile, jar):
	# Stores runtime options
	conf = {"ref":ref, "gatk":None, "picard":None}
	if cpu > cpu_count():
		cpu = cpu_count()
	conf["cpu"] = cpu
	conf["jar"] = jar
	if not os.path.isfile(conf["ref"]):
		print("\n\t[Error] Genome fasta not found. Exiting.\n")
		quit()
	if jar == True:
		print("\n\tReading config file...")
		with open(infile, "r") as f:
			for line in f:
				s = line.split("=")
				target = s[0].strip()
				val = s[1].strip()
				if val:
					if target == "GATK_jar":
						conf["gatk"] = val
					elif target == "Picard_jar":
						conf["picard"] = val
			# Check for errors
			if not conf["gatk"] or not os.path.isfile(conf["gatk"]):
				print("\n\t[Error] Please include path to GATK jar in config file. Exiting.\n")
				quit()
			if not conf["picard"] or not os.path.isfile(conf["picard"]):
				print("\n\t[Warning] Path to Picard.jar not found. Picard is required to create a new fasta dict.")
				proceed = input("\tProceed? (Enter 'y' for yes or any other key for no): ")
				if proceed.lower() != "y":
					print("\tExiting.\n")
					quit()
	return conf

def main():
	starttime = datetime.now()
	jar = False
	parser = argparse.ArgumentParser(description = "This script will call MuTect2 on a given \
list of input files. Be sure that pysam is installed and that Samtools is in your PATH.")
	parser.add_argument("-t", type = int, default = 1,
help = "The number of threads to use (i.e. number of parallel mutect instances).")
	parser.add_argument("-i", 
help = "Path to space/tab/comma seperated text file of input files (format: ID Normal A B)")
	parser.add_argument("-r", help = "Path to reference genome.")
	parser.add_argument("-o", help = "Path to output directory.")
	parser.add_argument("-j", help = "Path to config file. Will use the jar files indicated \
(calls binaries in environment by default).")
	args = parser.parse_args()
	if not args.i or not args.o or not args.r:
		print("\n\t[Error] Please specify input file, reference genome, and output directory. Exiting.\n")
		quit()
	if args.o[-1] != "/":
		args.o += "/"
	if args.j:
		jar = True
	conf = getConf(args.t, args.r, args.j, jar)
	checkReferences(conf)
	log, done, outfiles = checkOutput(args.o)
	conf["log"] = log
	files = getManifest(done, args.i)
	pool = Pool(processes = args.t)
	func = partial(submitFiles, conf, outfiles, args.o)
	l = len(files)
	# Call mutect
	print(("\n\tCalling mutect2 with {} threads....").format(args.t))
	for x in pool.imap_unordered(func, files):
		l -= 1
		print(("\n\t{} has finished. {} samples remaining").format(x, l))
	pool.close()
	pool.join()
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
