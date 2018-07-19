'''This script will mutect2 in parallel in tumor only mode and assemble a panel of normals'''

import os
from argparse import ArgumentParser
from datetime import datetime
from subprocess import Popen
from shlex import split
from runMutect import callMutect, getStatus
from bamUtil import checkRG

def makePON(infiles, outfile, gatk):
	# Calls mutect to create a new panel of normals
	first = True
	if gatk:
		# Format command for calling gatk jar
		cmd = ("java -jar {} CreateSomaticPanelOfNormals -O {}").format(gatk, outfile)
	else:
		cmd = ("gatk CreateSomaticPanelOfNormals -O {}").format(outfile)
	with open(infiles, "r") as f:
		for line in f:
			if first == False and line[0] != "#":
				splt = line.split("\t")
				if len(splt) >= 2:
					cmd += (" -vcfs {}").format(splt[1].strip())
			else:
				first = False
	log = outfile[:outfile.find(".")] + ".stdout"
	with open(log, "w") as l:
		try:
			mp = Popen(split(cmd), stdout = l, stderr = l)
			mp.communicate()
		except:
			print("\n\t[Error] Could not generate panel of normals. Exiting.")
			quit()

	return getStatus(log)

def submitNormal(conf):
	# Builds mutect command
	if "picard" in conf.keys():
		tumorname, bam = checkRG(conf["normal"], conf["sample"], conf["picard"])
	else:
		tumorname, bam = checkRG(conf["normal"], conf["sample"])
	if not tumorname or not bam:
		print("\t[Error] Failed adding read groups. Exiting")
		quit()
	# Assemble command
	outfile = conf["outpath"] + conf["sample"] + ".vcf"
	if "gatk" in conf.keys():
		# Format command for calling gatk jar
		cmd = ("java -jar {} Mutect2 --disable-read-filter \
MateOnSameContigOrNoMappedMateReadFilter -R {} ").format(conf["gatk"], conf["reference"])
	else:
		# Format command for colling gatk from path
		cmd = ("gatk Mutect2 --disable-read-filter \
MateOnSameContigOrNoMappedMateReadFilter -R {} ").format(conf["reference"])
	cmd += ("--tumor-sample {} -I {} --output {}").format(tumorname, bam, outfile)
	if "bed" in conf.keys():
		cmd += (" -L {}").format(conf["bed"])
	if "germline" in conf.keys():
		cmd += (" --germline-resource {} --af-of-alleles-not-in-resource {}").format(conf["germline"], conf["af"])
	# Call mutect for control and tumor
	res = callMutect(cmd, conf["sample"], outfile)
	if res:
		print("\tFinished running Mutect2.")
		with open(conf["log"], "a") as log:
			log.write(("{}\t{}\n").format(conf["sample"], outfile))
		return True
	else:
		print("\t[Error] Failed running Mutect2. Exiting")
		quit()

def configEntry(conf, arg, key):
	# Returns dict with updated arg entry
	if not arg:
		print(("\n\t[Error] Please specify {}. Exiting.\n").format(arg))
		quit()
	else:
		conf[key] = arg
	return conf	

def getConfig(args):
	# Returns arguments as dict
	conf = {}
	conf["log"] = args.l
	if args.o[-1] != "/":
		args.o += "/"
	conf = configEntry(conf, args.s, "sample")
	conf = configEntry(conf, args.c, "normal")
	conf = configEntry(conf, args.r, "reference")
	conf = configEntry(conf, args.o, "outpath")
	if args.bed:
		conf["bed"] = args.bed
	if args.gatk:
		conf["gatk"] = args.gatk
	if args.picard:
		conf["picard"] = args.picard
	if args.g:
		if not args.af:
			print("\n\t[Error] Please supply an allele frequency when using a germline estimate. Exiting.\n")
			quit()
		else:
			conf["germline"] = args.g
			conf["af"] = args.af
	if args.e:
		conf["contaminant"] = args.e	
	return conf

def main():
	starttime = datetime.now()
	parser = ArgumentParser(description = "This script will mutect2 in parallel in tumor \
only mode and assemble a panel of normals. Be sure that pysam is installed and that bcftools is in your PATH.")
	parser.add_argument("--pon", default = False, action = "store_true", help = "Generate new \
panel of normals from log file (requires -l (output from tumor only mode) and -o (output PON file) flags only).")
	parser.add_argument("-s", help = "Sample name (required).")
	parser.add_argument("-l", help = "Path to log file (required; output files are recorded here).")
	parser.add_argument("-c", help = "Path to normal/control bam (required).")
	parser.add_argument("-r", help = "Path to reference genome (required).")
	parser.add_argument("-o", help = "Path to output directory (required).")
	parser.add_argument("--bed", help = "Path to bed annotation.")
	parser.add_argument("--gatk", help = "Path to gatk jar (if using).")
	parser.add_argument("--picard", help = "Path to picard jar (if using).")
	parser.add_argument("-g", help = "Path to germline resource.")
	parser.add_argument("--af", help = "Estimated allele frequency (required if using a germline resource).")
	parser.add_argument("-e", help = "Path to contmination estimate vcf.")
	args = parser.parse_args()
	if args.pon == True:
		print("\n\tGenerating panel of normals...")
		status = makePON(args.l, args.o, args.gatk)
		if status == False:
			print("\n\t[Error] Could not generate panel of normals. Exiting.")
	else:
		conf = getConfig(args)
		# Call mutect
		print(("\n\tCalling Mutect2 in tumor-only mode on {}....").format(conf["sample"]))
		status = submitNormal(conf)
	if status == True:
		print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
