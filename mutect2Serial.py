'''This script will call MuTect2 on a given list of input files'''

import argparse
import os
import pysam
from datetime import datetime
from collections import OrderedDict
from subprocess import Popen
from shlex import split

def getTumorSample(bam):
	# Extracts tumor sample name from bam file
	rg = pysam.view("-H", bam)
	rg = rg[rg.find("@RG"):rg.find("@PG")]
	name = rg[rg.find("SM:"):]
	return name[name.find(":")+1:name.find("\t")]

def callMutect(cmd, path, n, t):
	# Calls Mutect with given root command and files
	nid = n[n.rfind("/")+1:].replace(".bam", "")
	tid = t[t.rfind("/")+1:].replace(".bam", "")
	outfile = ("{}-{}_{}.vcf").format(path, nid, tid)
	tumorname = getTumorSample(t)
	cmd += ("--tumor-sample {} -I {} -I {} --output {}").format(tumorname, t, n, outfile)
	print(("\tComparing {} and {}...").format(nid, tid))
	with open(os.devnull, "w") as dn:
		try:
			mt = Popen(split(cmd), stdout = dn)	
			mt.communicate()
			print("\tFinished.")
		except:
			print(("\t[Error] Could not call MuTect on {} and {}.").format(nid, tid))

def samIndex(bam):
	# Call samtools to index sam/bam file
	if not os.path.isfile(bam + ".bai"):
		print(("\tGenerating sam index for {}...").format(bam))
		pysam.index(bam)

def submitFiles(conf, files, outdir, jar):
	# Calls MuTect2 serially over input files
	print("\n\tCalling MuTect2 on input files...")
	for i in files.keys():
		print(("\n\tProcessing {}...").format(i))
		# Assemble command
		if jar == False:
			# Format command for colling gatk from path
			cmd = "gatk Mutect2 "
		else:
			# Format command for calling gatk jar
			cmd = ("java -jar {} Mutect2 ").format(conf["gatk"])
		cmd += ("-R {} ").format(conf["ref"])
		#cmd += ("-R {} --germline-resource {} -pon {} ").format(conf["ref"], conf["cosmic"], conf["dbsnp"])
		# Call for each combination of files
		samIndex(files[i][0])
		for j in files[i][1:]:
			samIndex(j)
			callMutect(cmd, outdir + i, files[i][0], j)
		# Record finished samples
		'''with open(conf["log"], "a") as l:
			l.write(i + "\n")'''
		print(("\tFinished running {}.").format(i))

def getManifest(done, infile):
	# Returns dict of input files
	files = OrderedDict()
	first = True
	print("\tReading input file...")
	with open(infile, "r") as f:
		for line in f:
			if first == True:
				# Determine delimiter
				for i in [" ", "\t", ","]:
					if i in line:
						delim = i
						break
				first = False
			s = line.strip().split(delim)
			if s[0] not in done:
				# {ID: [Normal, A, B]}
				files[s[0]] = s[1:]
	print(("\tFound entries for {} new samples.").format(len(files.keys())))
	return files

def checkOutput(outdir):
	# Checks for output log file and reads if present
	done = []
	log = outdir + "mutectLog.txt"
	print("\tChecking for previous output...")
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	if os.path.isfile(log):
		with open(log, "r") as f:
			for line in f:
				done.append(line.strip())
		print(("\tFound {} completed samples.").format(len(done)))
	else:
		with open(log, "w") as f:
			# Initialize log file
			pass
	return log, done

def checkReferences(conf, jar):
	# Ensures fasta and vcf index and dict files are present
	ref = conf["ref"]
	if not os.path.isfile(ref + ".fai"):
		# Call samtools to make index file
		print("\tGenerating fasta index...")
		fai = Popen(split(("samtools faidx {}").format(ref)))
		fai.communicate()
	fdict = ref.replace(".fa", ".dict")
	with open(os.devnull, "w") as dn:
		if not os.path.isfile(fdict):
			# Call CreateSequenceDictionary
			print("\tGenerating fasta dictionary...\n")
			if jar == True:
				cmd = ("java -jar {} ").format(conf["picard"])
			else:
				cmd = "picard "
			fd = Popen(split(("{} CreateSequenceDictionary R= {} O= {}").format(cmd, ref, fdict)), stdout=dn)
			fd.communicate()
		# Check for vcf indexes
		if jar == True:
			cmd = cmd = ("java -jar {} ").format(conf["gatk"])
		else:
			cmd = "gatk "
		cmd += "IndexFeatureFile -F "
		for i in [conf["cosmic"], conf["dbsnp"]]:
			if not os.path.isfile(i + ".idx"):
				print(("\tGenerating vcf index for {}...\n").format(i))
				vcfi = Popen(split(("{} {}").format(cmd, i)), stdout=dn)
				vcfi.communicate()

def getConf(infile, jar):
	# Stores runtime options
	conf = {"cosmic":None, "dbsnp":None, "ref":None, "gatk":None, "picard":None}
	print("\n\tReading config file...")
	with open(infile, "r") as f:
		for line in f:
			s = line.split("=")
			target = s[0].strip()
			val = s[1].strip()
			if val:
				if target == "COSMIC":
					conf["cosmic"] = val
				elif target == "dbSNP":
					conf["dbsnp"] = val
				elif target == "reference_genome":
					conf["ref"] = val
				elif target == "GATK_jar":
					conf["gatk"] = val
				elif target == "Picard_jar":
					conf["picard"] = val
		# Check for errors
		if not conf["cosmic"] or not conf["dbsnp"]:
			print("\t[Error] Please add paths to COSMIC and dbSNP vcf files to config file. Exiting.\n")
			quit()
		if jar == True:
			if not conf["gatk"]:
				print("\t[Error] Please include path to GATK jar in config file. Exiting.\n")
				quit()
			if not conf["picard"]:
				print("\t[Warning] Path to Picard.jar not found. Picard is required to create a new fasta dict.")
				proceed = input("\tProceed? (Enter 'y' for yes or any other key for no): ")
				if proceed.lower() != "y":
					print("\tExiting.\n")
					quit()
		if not conf["ref"]:
			print("\t[Error] Please include path to reference genome in config file. Exiting.\n")
			quit()
	return conf

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser(description = "This script will call MuTect2 on a given \
list of input files. Be sure that pysam is installed and that Samtools is in your PATH.")
	parser.add_argument("-i", 
help = "Path to space/tab/comma seperated text file of input files (format: ID Normal A B)")
	parser.add_argument("-c", help = "Path to config file.")
	parser.add_argument("-o", help = "Path to output directory.")
	parser.add_argument("--jar", action = "store_true", default = False,
help = "Use if not loading a module/environment (path to gatk and picard jars must be in config file.)")
	args = parser.parse_args()
	if not args.i or not args.o or not args.c:
		print("\n\t[Error] Please specify input file, config file, and output directory. Exiting.\n")
		quit()
	if args.o[-1] != "/":
		args.o += "/"
	conf = getConf(args.c, args.jar)
	checkReferences(conf, args.jar)
	log, done = checkOutput(args.o)
	conf["log"] = log
	files = getManifest(done, args.i)
	submitFiles(conf, files, args.o, args.jar)
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
