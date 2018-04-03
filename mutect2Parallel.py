'''This script will call MuTect2 on a given list of input files'''

import argparse
import os
from datetime import datetime
from functools import partial
from multiprocessing import Pool, cpu_count
from subprocess import Popen
from shlex import split
from bamUtil import *

def intersect(outfile, cmd):
	# Calls bcftools to get intersecting rows and summarizes output
	try:
		bcf = Popen(split(cmd))
		bcf.communicate()
	except:
		print(("\t[Error] Could not call bcftools with {}.").format(cmd))
		return 0
	# Initialize output csv
	#with open(outfile, "w") as output:
	#	output.write("SampleA,SampleB,#PrivateA,#PrivateB,#Common,Similarity\n")
	return 1

def filterCalls(cmd, vcf):
	# Calls gatk to filter mutect calls
	outfile = vcf[:vcf.find(".")] + ".filtered.vcf"
	outfile = outfile.replace("Mutect/", "Filtered/")
	log = outfile.replace("vcf", "stdout")
	cmd += ("-V {} -O {}").format(vcf, outfile)
	with open(log, "w") as l:
		try:
			fmc = Popen(split(cmd), stdout = l, stderr = l)
			fmc.communicate()
		except:
			print(("\t[Error] Could not call FilterMutectCalls on {}.").format(vcf))
			return ""
	# bgzip and index output
	idxfile = tabixIndex(vcf)
	return idxfile

def compareVCFs(conf, outfile, vcfs):
	# Calls gatk and pyvcf to filter and compare mutect output
	ret = False
	filtered = []
	if conf["jar"] == False:
		# Format command for colling gatk from path
		cmd = "gatk FilterMutectCalls "
	else:
		# Format command for calling gatk jar
		cmd = ("java -jar {} FilterMutectCalls ").format(conf["gatk"])
	for i in vcfs:
		# Filter each vcf
		x = filterCalls(cmd, i)
		if x:
			filtered.append(x)
	if len(filtered) == 2:
		done = 0
		# Call bftools on all results
		cmd = ("bcftools isec {} {} -p {}").format(filtered[0], filtered[1])
		done += intersect(outfile, cmd)
		# Call bftools on passes
		done += intersect(outfile.replace(".csv", "_PASS.csv"), cmd + " -f .,PASS")
		if done == 2:
			ret = True
	return ret

def callMutect(cmd, picard, path, n, t):
	# Calls Mutect with given root command and files
	nid = n[n.rfind("/")+1:].replace(".sorted", "").replace(".bam", "")
	tid = t[t.rfind("/")+1:].replace(".sorted", "").replace(".bam", "")
	outfile = ("{}-{}.vcf").format(path, tid)
	tumorname, t = checkRG(t, tid, picard)
	if not tumorname:
		return ""
	# Build command and call mutect
	cmd += ("--tumor-sample {} -I {} -I {} --output {}").format(tumorname, t, n, outfile)
	print(("\tComparing {} and {}...").format(nid, tid))
	# Make log file
	log = outfile.replace("vcf", "stdout")
	with open(log, "w") as dn:
		try:
			mt = Popen(split(cmd), stdout = dn, stderr = dn)	
			mt.communicate()
		except:
			print(("\t[Error] Could not call MuTect on {} and {}.").format(nid, tid))
			return ""
	with open(log, "r") as dn:
		status = False
		# Make sure mutect completed successfully
		for line in dn:
			if "Tool returned:" in line:
				status = True
			elif status == True:
				# Get exit status
				if "SUCCESS" in line:
					print(("\tFinished comparing {} and {}...").format(nid, tid))
					return outfile
				else:
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
		# Get list of bam files that need to be run
		if sample[4] == 1:
			s = [sample[3]]
		elif sample[4] == 2:
			s = [sample[2]]
		else:
			s = [sample[2], sample[3]]
		# Call for each remaining combination of control-tumor
		for j in s:
			res = callMutect(cmd, conf["picard"], outdir + "Mutect/" + sample[0], control, j)
			if res:
				# Record finished samples
				vcfs.append(res)
				with open(conf["log"], "a") as l:
					l.write(("{}\tMutect\t{}\n").format(res[res.rfind("/")+1:res.rfind(".")], res))
	if sample[4] == 3 or len(vcfs) == 2:
		# Compare output
		print(("\n\tComparing VCFs from {}...").format(sample[0]))
		status = compareVCFs(conf, outdir + sample[0] + ".csv", vcfs)
		if status == True:
			# Record finished samples
			with open(conf["log"], "a") as l:
				l.write(("{}\tcomparison\n").format(sample[0]))
	if len(vcfs) < 2:
		return ("Could not complete {}").format(sample[0])
	else:
		return ("{}  has finished").format(sample[0])

#-------------------------------------------Inputs----------------------------

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
					if s[2][s[2].rfind("/")+1:s[2].find(".")] in done[s[0]]:
						a = 1
					if s[3][s[3].rfind("/")+1:s[3].find(".")] in done[s[0]]:
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
	if not os.path.isdir(outdir + "Mutect/"):
		os.mkdir(outdir + "Mutect/")
	if not os.path.isdir(outdir + "Filtered/"):
		os.mkdir(outdir + "Filtered/")
	if os.path.isfile(log):
		with open(log, "r") as f:
			for line in f:
				if first == False and line.strip():
					line = line.strip().split("\t")
					splt = line[0].split("-")
					if len(line) == 3 and len(splt) == 2:
						if splt[0] in done.keys():
							# Record vcf location
							outfiles[splt[0]].append(line[2])
							if line[1] == "comparison":
								# Record last stage completed (comparison/mutect)
								done[splt[0]] = [done[splt[0]][0], splt[1], line[1]]
							else:
								done[splt[0]] = [done[splt[0]][0], splt[1], done[splt[0]][1]]
						else:
							# {sample#: [filename1, filename2, mostRecent]}
							done[splt[0]] = [splt[1], line[1]]
							# {sample#: [vcf1, vcf2]}
							outfiles[splt[0]] = [line[2]]
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
	conf = {"ref":None, "gatk":None, "picard":None}
	if cpu > cpu_count():
		cpu = cpu_count()
	conf["cpu"] = cpu
	conf["jar"] = jar
	if ref:
		conf["ref"] = ref
	if not ref or jar == True:
		print("\n\tReading config file...")
		with open(infile, "r") as f:
			for line in f:
				s = line.split("=")
				target = s[0].strip()
				val = s[1].strip()
				if val:
					if target == "reference_genome":
						conf["ref"] = val
					elif target == "GATK_jar":
						conf["gatk"] = val
					elif target == "Picard_jar":
						conf["picard"] = val
			# Check for errors
			if jar == True:
				if not conf["gatk"] or not os.path.isfile(conf["gatk"]):
					print("\n\t[Error] Please include path to GATK jar in config file. Exiting.\n")
					quit()
				if not conf["picard"] or not os.path.isfile(conf["picard"]):
					print("\n\t[Error] Path to Picard.jar not found. Exiting.\n")
					quit()
	if not os.path.isfile(conf["ref"]):
		print("\n\t[Error] Genome fasta not found. Exiting.\n")
		quit()
	return conf

#-----------------------------------------------------------------------------

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser(description = "This script will call MuTect2 on a given \
list of input files. Be sure that pysam is installed and that bcftools is in your PATH.")
	parser.add_argument("--jar", action = "store_true", default = False,
help = "Calls java jars indicated in the config file (Calls binaries in the environment by default).")
	parser.add_argument("-t", type = int, default = 1,
help = "The number of threads to use (i.e. number of parallel mutect instances).")
	parser.add_argument("-i", 
help = "Path to space/tab/comma seperated text file of input files (format: ID Normal A B)")
	parser.add_argument("-r", help = "Path to reference genome.")
	parser.add_argument("-o", help = "Path to output directory.")
	parser.add_argument("-c", 
help = "Path to config file containing reference genome, java jars (if using), and mutect options).")
	args = parser.parse_args()
	if not args.i or not args.o:
		print("\n\t[Error] Please specify input file and output directory. Exiting.\n")
		quit()
	if args.o[-1] != "/":
		args.o += "/"
	conf = getConf(args.t, args.r, args.c, args.jar)
	checkReferences(conf)
	log, done, outfiles = checkOutput(args.o)
	conf["log"] = log
	files = getManifest(done, args.i)
	pool = Pool(processes = args.t)
	func = partial(submitFiles, conf, outfiles, args.o)
	l = len(files)
	t = "threads"
	if args.t > 1:
		t += "s"
	# Call mutect+
	print(("\n\tCalling mutect2 with {} {}....").format(args.t, t))
	for x in pool.imap_unordered(func, files):
		l -= 1
		print(("\n\t{}. {} samples remaining").format(x, l))
	pool.close()
	pool.join()
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
