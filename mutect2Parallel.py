'''This script will call MuTect2 on a given list of input files'''

import os
from argparse import ArgumentParser
from datetime import datetime
from subprocess import Popen
from shlex import split
from bamUtil import getFastaIndex

def submitJobs(scripts, batch):
	# Determines grid type and submits jobs
	cmd = ""
	for line in batch:
		if "#SBATCH" in line:
			cmd = "sbatch "
			break
		elif "#PBS" in line:
			cmd = "qsub "
			break
	if len(cmd) <= 3:
		print("\t[Error] Cannot determine grid type. Exiting.\n")
		quit()
	for i in scripts:
		try:
			# Submit each batch script
			Popen(split(cmd + i + " \n"))
		except:
			print(("\t[Error] Could not submit {}").format(i))
	print()
	return True

def getCommand(conf):
	# Returns base python call for all files
	if conf["newpon"] == False:
		# Format tumor-normal run
		cmd = ("python runPair.py -r {} ").format(conf["ref"])
		if conf["nofilter"] == True:
			cmd += "--nofilter "
		elif conf["filter"] == True:
			cmd += "--filter "
		if conf["bamout"] == True:
			cmd += "--bamout "
		if "pon" in conf.keys():
			cmd += ("-p {} ").format(conf["pon"])
		if "germline" in conf.keys():
			cmd += ("-g {} ").format(conf["germline"])
		if "contaminant" in conf.keys() and "af" in conf.keys():
			cmd += ("--af {} -e {} ").format(conf["af"], conf["contaminant"])
		if "mo" in conf.keys():
			cmd += ('--mo "{} "').format(conf["mo"])
		if "fmo" in conf.keys():
			cmd += ('--fmo "{} "').format(conf["fmo"])
	else:
		# Format for tumor-only mode
		cmd = ("python getPON.py -l {} -r {} ").format(conf["outpath"] + "normalsLog.txt", conf["ref"])
		if not os.path.isfile(conf["outpath"] + "normalsLog.txt"):
			with open(conf["outpath"] + "normalsLog.txt", "w") as f:
				# initilize lof file
				f.write("Sample\tVCF\n")
	for i in ["bed", "gatk", "picard"]:
		if i in conf.keys() and conf[i] != None:
			cmd += ("--{} {} ").format(i, conf[i])
	return cmd

def getBatchScripts(outdir, conf, batch, files):
	# Generates new batch script for each set of samples
	cmd = getCommand(conf)
	scripts = []
	for i in files.keys():
		print(("\tWriting batch script for {}...").format(i))
		outfile = outdir + i + ".sh"
		with open(outfile, "w") as output:
			for line in batch:
				if "--job-name=" in line:
					output.write(("{}_{}\n").format(line.strip(), i))
				else:
					output.write(line)
			if conf["newpon"] == False:
				outpath = conf["outpath"] + i + "/"
				c = cmd + ("-s {} -c {} -x {} -y {} -o {}").format(i, 
					files[i][0], files[i][1], files[i][2], outpath)
			else:
				c = cmd + ("-s {} -c {} -o {}").format(i, files[i], conf["outpath"])
			output.write(c + "\n")
		scripts.append(outfile)
	return scripts

def getManifest(infile, pon):
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
				if pon == False:
					files[s[0]] = s[1:]
				else:
					# Only save normal file path
					files[s[0]] = s[1]
	for i in files.keys():
		if pon == False:
			for j in files[i]:
				if not os.path.isfile(j):
					print(("\n\t[Error] Input file {} not found. Exiting.\n").format(j))
					quit()
		else:
			if not os.path.isfile(files[i]):
				print(("\n\t[Error] Input file {} not found. Exiting.\n").format(files[i]))
				quit()
	print(("\tFound entries for {} samples.\n").format(len(files)))
	return files

#-----------------------------------------------------------------------------

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
			if "picard" in conf.keys():
				cmd = ("java -jar {} ").format(conf["picard"])
			else:
				cmd = "picard "
			fd = Popen(split(("{} CreateSequenceDictionary R= {} O= {}").format(cmd, ref, fdict)), stdout=dn, stderr=dn)
			fd.communicate()
		for i in ["pon", "contaminant"]:
			if i in conf.keys():
				if not os.path.isfile(conf[i] + ".tbi"):
					if "gatk" in conf.keys():
						cmd = ("java -jar {} IndexFeatureFile -F {}").format(conf["gatk"], conf[i])
					else:
						cmd = ("gatk IndexFeatureFile -F {}").format(conf[i])
					fd = Popen(split(cmd), stdout=dn, stderr=dn)
					fd.communicate()

def checkFile(infile, err):
	# Exits if infile is not found
	if not os.path.isfile(infile):
		print(("\n\t[Error] {} not found. Exiting.\n").format(err))
		quit()	

def getOptions(conf, line):
	# Returns config dict with updated options
	val = None
	s = line.split("=")
	if len(s) == 2:
		target = s[0].strip()
		val = s[1].strip()
	if val:
		if target == "reference_genome":
			conf["ref"] = val
		elif target == "bed_annotation":
			conf["bed"] = val
			checkFile(conf["bed"], "Bed file")
		elif target == "output_directory":
			if val[-1] != "/":
				val += "/"
			conf["outpath"] = val
		elif target == "GATK_jar":
			conf["gatk"] = val
			checkFile(conf["gatk"], "GATK jar")
		elif target == "Picard_jar":
			conf["picard"] = val
			checkFile(conf["picard"], "Picard jar")	
		elif target == "normal_panel":
			conf["pon"] = val
			checkFile(conf["pon"], "Panel of normals file")
		elif target == "germline_resource":
			conf["germline"] = val
			checkFile(conf["germline"], "Germline resource file")
		elif target == "allele_frequency":
			conf["af"] = val
		elif target == "contaminant_estimate":
			conf["contaminant"] = val
			checkFile(conf["contaminant"], "Contaminant estimate vcf")
		elif target == "mutect_options":
			# Extract directly from line in case options have an equals sign
			conf["mo"] = line[line.find("=")+1:]
			conf["mo"] = conf["mo"].strip()
		elif target == "filter_mutect_options":
			conf["fmo"] = line[line.find("=")+1:]
			conf["fmo"] = conf["mo"].strip()	
	return conf

def getConf(infile):
	# Stores runtime options
	batch = []
	opt = True
	store = False
	conf = {"ref":None}
	print("\n\tReading config file...")
	with open(infile, "r") as f:
		for line in f:
			if "Sample Batch Script" in line:
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
	# Check for critical errors
	if not os.path.isfile(conf["ref"]):
		print("\n\t[Error] Genome fasta not found. Exiting.\n")
		quit()
	if not os.path.isdir(conf["outpath"]):
		os.mkdir(conf["outpath"])
	if "germline" in conf.keys() and "af" not in conf.keys():
		print("\n\t[Error] Please supply an allele frequency when using a germline estimate. Exiting.\n")
		quit()	
	return conf, batch

def main():
	starttime = datetime.now()
	parser = ArgumentParser(description = "This script will call MuTect2 on a given \
list of input files. Be sure that pysam is installed and that bcftools is in your PATH.")
	parser.add_argument("--submit", action = "store_true", default = False,
help = "Submit batch files to SLURM/Torque grid for execution.")
	parser.add_argument("--bamout", action = "store_true", default = False,
help = "Indicates that mutect should also generate bam output files.")
	parser.add_argument("--newPON", action = "store_true", default = False,
help = "Creates batch scripts for running mutect on normals and creating a panel of normals.")
	parser.add_argument("--nofilter", action = "store_true", default = False,
help = "Skips filtering of mutect output.")
	parser.add_argument("--filter", action = "store_true", default = False,
help = "Begins pipeline with filtering of previous mutect output")
	parser.add_argument("-i", 
help = "Path to space/tab/comma seperated text file of input files (format: ID Normal A B)")
	parser.add_argument("-c", 
help = "Path to config file containing reference genome, java jars (if using), and mutect options.")
	parser.add_argument("-o", default = "",
help = "Path to batch script output directory (leave blank for current directory).")
	args = parser.parse_args()
	if not args.i and args.c:
		print("\n\t[Error] Please specify input file and config file. Exiting.\n")
		quit()
	if args.nofilter == True and args.filter == True:
		print("\n\t[Error] Please specify only one of filter/nofilter. Exiting.\n")
		quit()
	if args.o and args.o[-1] != "/":
		args.o += "/"
	conf, batch = getConf(args.c)
	conf["bamout"] = args.bamout
	conf["newpon"] = args.newPON
	conf["nofilter"] = args.nofilter
	conf["filter"] = args.filter
	checkReferences(conf)
	files = getManifest(args.i, conf["newpon"])
	scripts = getBatchScripts(args.o, conf, batch, files)
	if args.submit == True:
		done = submitJobs(scripts, batch)
	else:
		done = True
	if done == True:
		print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
