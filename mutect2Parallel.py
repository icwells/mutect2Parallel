'''This script will call MuTect2 on a given list of input files'''

import os
from argparse import ArgumentParser
from datetime import datetime
from subprocess import Popen
from shlex import split
from bamUtil import getFastaIndex

def getCommand(conf):
	# Returns base python call for all files
	cmd = ("python runPair.py -r {} ").format(conf["ref"])
	for i in ["bed", "gatk", "picard"]:
		if i in conf.keys() and conf[i] != None:
			cmd += ("-{} {} ").format(i, conf[i])
	return cmd

def getBatchScripts(outdir, path, conf, batch, files):
	# Generates new batch script for each set of samples
	count = 1
	cmd = getCommand(conf)
	for i in files.keys():
		print(("\tWriting batch script for {}...").format(i))
		with open(outdir + i + ".sh", "w") as output:
			for line in batch:
				if "--job-name=" in line:
					output.write(("{}_{}\n").format(line.strip(), count))
				else:
					output.write(line)
			outpath = path + i + "/"
			c = cmd + ("-s {} -c {} -x {} -y {} -o {}").format(i, 
					files[i][0], files[i][1], files[i][2], outpath)
			output.write(c + "\n")
		count += 1

def getManifest(infile):
	# Returns dict of input files
	files = {}
	first = True
	print("\tReading input file...")
	with open(infile, "r") as f:
		for line in f:
			if line[0] != "#":
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
	print(("\tFound entries for {} samples.\n").format(len(files)))
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

def getOptions(conf, line):
	# Returns config dict with updated options
	s = line.split("=")
	target = s[0].strip()
	val = s[1].strip()
	if val:
		if target == "reference_genome":
			conf["ref"] = val
		elif target == "bed_annotation":
			conf["bed"] = val
		elif target == "GATK_jar":
			conf["gatk"] = val
		elif target == "Picard_jar":
			conf["picard"] = val
	return conf

def getConf(infile):
	# Stores runtime options
	batch = []
	opt = True
	store = False
	conf = {"ref":None}
	#conf = {"ref":None, "gatk":None, "picard":None, "bed":None}
	print("\n\tReading config file...")
	with open(infile, "r") as f:
		for line in f:
			if line.count("=") >= 10:
				opt = False
			elif opt == True and line.strip():
				# Store options
					conf = getOptions(conf, line)
			else:
				if "#!/" in line:
					store = True
				if store == True:
					# Store batch script
					batch.append(line)
	# Check for errors
	if not os.path.isfile(conf["ref"]):
		print("\n\t[Error] Genome fasta not found. Exiting.\n")
		quit()
	for i in ["gatk", "picard", "bed"]:
		if i in conf.keys():
			if not os.path.isfile(conf[i]):
				if i == "bed":
					print("\n\t[Error] Bed annotation file not found. Exiting.\n")
				else:
					print(("\n\t[Error] {} jar not found. Exiting.\n").format(i.upper))
				quit()
	return conf, batch

#-----------------------------------------------------------------------------

def main():
	starttime = datetime.now()
	parser = ArgumentParser(description = "This script will call MuTect2 on a given \
list of input files. Be sure that pysam is installed and that bcftools is in your PATH.")
	parser.add_argument("-i", 
help = "Path to space/tab/comma seperated text file of input files (format: ID Normal A B)")
	parser.add_argument("-c", 
help = "Path to config file containing reference genome, java jars (if using), and mutect options).")
	parser.add_argument("-p", help = "Path to mutect output directory.")
	parser.add_argument("-o", default = "",
help = "Path to batch script output directory (leave blank for current directory).")
	args = parser.parse_args()
	if not args.i or not args.p:
		print("\n\t[Error] Please specify input file and output directory. Exiting.\n")
		quit()
	if args.o and args.o[-1] != "/":
		args.o += "/"
	if args.p[-1] != "/":
		args.p += "/"
	conf, batch = getConf(args.c)
	checkReferences(conf)
	files = getManifest(args.i)
	getBatchScripts(args.o, args.p, conf, batch, files)
	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
