'''This script will call MuTect2 on a given list of input files'''

import os
from argparse import ArgumentParser
from datetime import datetime
from subprocess import Popen
from shlex import split
from bamUtil import getFastaIndex

def getBatchScripts(conf, files):
	# Generates new batch script for each set of samples

def getManifest(infile):
	# Returns dict of input files
	files = {}
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
			files[s[0]] = s[1:]
	for i in files.keys():
		for j in files[i]:
			if not os.path.isfile(j):
				print(("\n\t[Error] Input file {} not found. Exiting.\n").format(j))
				quit()
	print(("\tFound entries for {} samples.").format(len(files)))
	return files

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

def getConf(cpu, infile, jar):
	# Stores runtime options
	conf = {"ref":None, "gatk":None, "picard":None, "mail":None, "bed":None "job":None}
	if cpu > cpu_count():
		cpu = cpu_count()
	conf["cpu"] = cpu
	conf["jar"] = jar
	if not ref or jar == True:
		print("\n\tReading config file...")
		with open(infile, "r") as f:
			for line in f:
				s = line.split("=")
				target = s[0].strip()
				val = s[1].strip()
				if val:
					if target == "e-mail":
						conf["mail"] = val
					elif target = "job_name":
						conf["job"] = val
					elif target == "reference_genome":
						conf["ref"] = val
					elif target == "bed_annotation":
						conf["bed"] = val
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
	if conf["bed"] != None:
		if not os.path.isfile(conf["bed"]):
			print("\n\t[Error] Bed annotation file not found. Exiting.\n")
			quit()
	return conf

#-----------------------------------------------------------------------------

def main():
	starttime = datetime.now()
	parser = argparse.ArgumentParser(description = "This script will call MuTect2 on a given \
list of input files. Be sure that pysam is installed and that bcftools is in your PATH.")
	parser.add_argument("--jar", action = "store_true", default = False,
help = "Calls java jars indicated in the config file (Calls binaries in the environment by default).")
	parser.add_argument("-i", 
help = "Path to space/tab/comma seperated text file of input files (format: ID Normal A B)")
	parser.add_argument("-o", help = "Path to output directory.")
	parser.add_argument("-c", 
help = "Path to config file containing reference genome, java jars (if using), and mutect options).")
	args = parser.parse_args()
	if not args.i or not args.o:
		print("\n\t[Error] Please specify input file and output directory. Exiting.\n")
		quit()
	if args.o[-1] != "/":
		args.o += "/"
	conf = getConf(args.t, args.c, args.jar)
	checkReferences(conf)
	files = getManifest(done, args.i, args.s)
	getBatchScripts(conf, files)
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
