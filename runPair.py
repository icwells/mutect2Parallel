'''This script will run a pair of input tumor bam files through mutect2 in parallel'''

import os
from argparse import ArgumentParser
from datetime import datetime
from functools import partial
from multiprocessing import Pool, cpu_count
from datetime import datetime
from subprocess import Popen
from shlex import split
from bamUtil import *

def getTotal(vcf):
	# Returns total number of content lines from vcf
	count = 0
	if os.path.isfile(vcf):
		with open(vcf, "r") as f:
			for line in f:
				if line[0] != "#":
					count += 1
	return count

def intersect(outpath, cmd, filtered):
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
	sa = filtered[0][filtered[0].rfind("/")+1:filtered[0].find(".")]
	sb = filtered[1][filtered[1].rfind("/")+1:filtered[1].find(".")]
	with open(outpath + ".csv", "w") as output:
		output.write("SampleA,SampleB,#PrivateA,#PrivateB,#Common,Similarity\n")
		output.write(("{},{},{},{},{},{}\n").format(sa, sb, a, b, c, sim))
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
			print(("\t[Error] Could not call FilterMutectCalls on {}").format(vcf))
			return ""
	return bgzip(outfile)

def compareVCFs(conf, outpath, vcfs):
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
		cmd = ("bcftools isec {} {}").format(filtered[0], filtered[1])
		cmd += " -p {}"
		done += intersect(outpath, cmd, filtered)
		# Call bftools on passes
		outpath += "_PASS"
		done += intersect(outpath, cmd + " -f .,PASS", filtered)
		if done == 2:
			ret = True
	return ret

#-------------------------------Mutect----------------------------------------

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
			print(("\t[Error] Could not call MuTect on {} and {}").format(nid, tid))
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
	if sample[0] in outfiles.keys():
		vcfs = outfiles[sample[0]]
	else:
		vcfs = []
	if sample[4] < 3:
		# Proceed if at least one combination has not been run
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
		status = compareVCFs(conf, outdir + sample[0], vcfs)
		if status == True:
			# Record finished samples
			with open(conf["log"], "a") as l:
				l.write(("{}\tcomparison\n").format(sample[0]))
	if len(vcfs) < 2:
		return ("Could not complete {}").format(sample[0])
	else:
		return ("{}  has finished").format(sample[0])

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

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser(description = "This script will call MuTect2 on a given \
list of input files. Be sure that pysam is installed and that bcftools is in your PATH.")
	parser.add_argument("-t", type = int, default = 2,
help = "The number of threads to use (default = 2: runs each member of a sample pair at once).")


	log, done, outfiles = checkOutput(args.o)
	conf["log"] = log
	pool = Pool(processes = args.t)
	func = partial(submitFiles, conf, outfiles, args.o)

	l = len(files)
	if l > 0:
		t = "thread"
		if args.t > 1:
			t += "s"
		# Call mutect+
		print(("\n\tCalling mutect2 with {} {}....").format(args.t, t))
		for x in pool.imap_unordered(func, files):
			l -= 1
			print(("\n\t{} has completed.").format(x))
		pool.close()
		pool.join()
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
