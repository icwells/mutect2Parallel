'''This script will filter mutect2 output files.'''

import os
from shlex import split
from commonUtil import *
from runPair import Sample, checkOutput, configEntry
from runMutect import appendLog, getSample

#-------------------------------Filtering-------------------------------------

def getTotal(vcf):
	# Returns total number of content lines from vcf
	count = 0
	if os.path.isfile(vcf):
		with open(vcf, "r") as f:
			for line in f:
				if line[0] != "#":
					count += 1
	return count

def intersect(outpath, cmd, vcfs):
	# Calls bcftools to get intersecting rows and summarizes output
	cmd = cmd.format(outpath)
	if "PASS" in outpath:
		typ = "pass"
	else:
		typ = "all"
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
	sa = vcfs[0][vcfs[0].rfind("/")+1:vcfs[0].find(".")]
	sb = vcfs[1][vcfs[1].rfind("/")+1:vcfs[1].find(".")]
	with open(outpath.replace("_PASS", "") + ".csv", "a") as output:
		output.write(("{},{},{},{},{},{},{}\n").format(typ,sa, sb, a, b, c, sim))
	return 1

def compareVCFs(conf, vcfs):
	# Calls gatk and pyvcf to filter and compare mutect output
	ret = False
	done = 0
	outpath = conf["outpath"] + conf["sample"]
	with open(outpath + "csv", "w") as output:
		# Initialize summary file and write header
		output.write("Type,SampleA,SampleB,#PrivateA,#PrivateB,#Common,Similarity\n")
	# Call bftools on all results
	cmd = ("bcftools isec {} {}").format(vcfs[0], vcfs[1])
	cmd += " -p {}"
	done += intersect(outpath, cmd, vcfs)
	# Call bftools on passes
	outpath += "_PASS"
	done += intersect(outpath, cmd + " -f .,PASS", vcfs)
	if done == 2:
		ret = True
	return ret

def filterCalls(conf, vcf):
	# Calls gatk to filter mutect calls
	outfile = vcf[:vcf.find(".")] + ".filtered.vcf"
	log = outfile.replace("vcf", "stdout")
	# Assemble command
	if "gatk" in conf.keys():
		cmd = ("java -jar {} FilterMutectCalls ").format(conf["gatk"])
	else:
		cmd = "gatk FilterMutectCalls "
	if "contaminant" in conf.keys():
		cmd += ("-contamination-table {} ").format(conf["contaminant"])
	if "fmo" in conf.keys():
		cmd += " " + conf["fmo"]
	cmd += ("-V {} -O {}").format(vcf, outfile)
	with open(log, "w") as l:
		try:
			fmc = Popen(split(cmd), stdout = l, stderr = l)
			fmc.communicate()
		except:
			print(("\t[Error] Could not call FilterMutectCalls on {}").format(vcf))
			return None
	return bgzip(outfile)

#-------------------------------Contamination---------------------------------

def estContamination(conf, s):
	# Calls gatk to get pileup summary and estimate contamination
	print(("\tEstimating contamination in {}...").format(s.ID))
	if "gatk" in conf.keys():
		# Format command for calling gatk jar
		pu = ("java -jar {} GetPileupSummaries -R {} ").format(conf["gatk"], conf["reference"])
		cc = ("java -jar {} CalculateContamination ").format(conf["gatk"])
	else:
		# Format command for colling gatk from path
		pu = ("gatk GetPileupSummaries -R {} ").format(conf["reference"])
		cc = "gatk CalculateContamination "
	# Get pileup summary
	pileup = s.Output.replace(".vcf", ".pileup.table")
	pu += ("-I {} -V {} -O {}").format(s.Input, conf["contaminant"], pileup)
	pu = getOpt(conf, pu)
	plog = pileup.replace("table", "stdout")
	with open(plog, "w") as l:
		try:
			spu = Popen(split(pu), stdout = l, stderr = l)
			spu.communicate()
		except:
			print(("\t[Error] Could not call GetPileupSummaries on {}").format(s.Output))
			return False
	if getStatus(plog) == False:
		return False
	# Get contamination estimate
	cest = pileup.replace("pileup", "contamination")
	clog = cest.replace("table", "stdout")
	cc += ("-I {} -O {}").format(pileup, cest)
	with open(clog, "w") as l:
		try:
			ccont = Popen(split(cc), stdout = l, stderr = l)
			ccont.communicate()
		except:
			print(("\t[Error] Could not call CalculateContamination on {}").format(s.Output))
			return False
	return getStatus(clog)

#-----------------------------------------------------------------------------

def filterSamples():
	# Filters and estimates contamination
	if s.Status == "mutect":
		if "contaminant" in conf.keys():
			status = estContamination(conf, s)
			if status == True:
				s.Status = "contamination-estimate:complete"
			else:
				s.Status = "contamination-estimate:failed"
		else:
			s.Status = "contamination-estimate:none"
		appendLog(conf, s)
	if "contamination-estimate" in s.Status and conf["nofilter"] == False:
		# Filter vcf
		filtered = filterCalls(conf, s.Output)
		if filtered:
			# Record filtered reads
			s.Output = filtered
			s.Status = "filtered"
		else:
			s.Status = "failed-filtering"
		appendLog(conf, s)

#-------------------------------IO--------------------------------------------

def getConfig(args):
	# Returns arguments as dict
	conf = {}
	if args.o[-1] != "/":
		args.o += "/"
	conf = configEntry(conf, args.s, "sample")
	conf = configEntry(conf, args.x, "tumor1")
	conf = configEntry(conf, args.r, "reference")
	conf = configEntry(conf, args.o, "outpath")
	if args.y:
		conf = configEntry(conf, args.y, "tumor2")
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
	if args.mo:
		conf["fmo"] = args.fmo
	return conf

def main():
	starttime = datetime.now()
	parser = ArgumentParser("This script will filter mutect2 output files.")
	parser.add_argument("-s", help = "Sample name (required).")
	parser.add_argument("-x", help = "Path to first vcf (required).")
	parser.add_argument("-y", help = "Path to second vcf (required).")
	parser.add_argument("-r", help = "Path to reference genome (required).")
	parser.add_argument("-o", help = "Path to output directory (required).")
	parser.add_argument("--gatk", help = "Path to gatk jar (if using).")
	parser.add_argument("--picard", help = "Path to picard jar (if using).")
	parser.add_argument("-g", help = "Path to germline resource.")
	parser.add_argument("--af", help = "Estimated allele frequency (required if using a germline resource).")
	parser.add_argument("-e", help = "Path to contmination estimate vcf.")
	parser.add_argument("--fmo", help = "Additional filter mutect options in quotes")
	args = parser.parse_args()
	conf = getConfig(args)
	log, samples = checkOutput(conf["outpath"])
	conf["log"] = log

	if len(filtered) == 2:
		# Compare output
		print(("\n\tComparing filterd VCFs from {}...").format(conf["sample"]))
		status = compareVCFs(conf, filtered)
		if status == True:
			# Record finished samples
			outfile = conf["outpath"] + conf["sample"] + ".csv"
			with open(conf["log"], "a") as l:
				l.write(("{}\tcompleted\t{}\n").format(conf["sample"], outfile))

	print(("\n\tFinished. Runtime: {}\n").format(datetime.now()-starttime))

if __name__ == "__main__":
	main()
