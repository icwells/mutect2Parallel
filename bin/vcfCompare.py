'''This script contains fucntions for filtering and comparing vcf files'''

import os
from bamUtils import *

def printError(msg):
	# Prints formatted error message
	print(("\n\t[Error] {}. Skipping.\n").format(msg))

#-------------------------------Comparison------------------------------------

def comparePair(outpath, vcfs):
	# Calls gatk and pyvcf to filter and compare given pair of
	for i in range(len(vcfs)):
		if vcfs[i] and os.path.isfile(vcfs[i]):
			# Make sure files are bgzipped
			vcfs[i] = tabix(vcfs[i])
		elif ".gz" not in vcfs[i] and os.path.isfile(vcfs[i] + ".gz"):
			vcfs[i] += ".gz"
		else:
			printError(("Cannot find {}").format(vcfs[i]))
			return None
	# Call bftools on all results
	ret = bcfIsec(outpath, vcfs)
	return ret

def compareVCFs(conf, log, name, samples):
	# Compares unfilted vs. passed results for each combination of pair of samples
	ret = False
	outpath = conf["outpath"] + name + "/"
	s1, s2 = list(samples.keys())
	aout = outpath + samples[s1].ID
	bout = outpath + samples[s2].ID
	# Append filtered sample name when submitting
	if getTotal(samples[s1].Output) > 0:
		a = comparePair(aout, [samples[s1].Output, samples[s2].Unfiltered])
	else:
		a = 0
	if getTotal(samples[s2].Output) > 0:
		b = comparePair(bout, [samples[s2].Output, samples[s1].Unfiltered])
	else:
		b = 0
	if a > 0 and b > 0:
		# Merge common variants and get total and similarity
		acom = tabix(aout + "/0002.vcf")
		bcom = tabix(bout + "/0002.vcf")
		common = bcfMerge(outpath, [acom, bcom])
		if common and os.path.isfile(common):
			c = getTotal(common)
			try:
				sim = c/(a+b+c)
			except ZeroDivisionError:
				sim = 0.0
	else:
		c = 0
		sim = 0.0
	with open(log, "a") as out:
		out.write(("{},{},{},{},{},{},{:.2%}\n").format(name, samples[s1].ID, samples[s2].ID, a, b, c, sim))
	ret = True
	return ret

#-------------------------------Filtering-------------------------------------

def snpsiftOpt(conf, typ):
	# Formats command with approriate parameters
	k = conf.keys()
	params = "("
	if typ == "a":
		params += "(FILTER != germline_risk) & "
		if "qual" in k:
			params += "(QUAL >= {}) & ".format(conf["qual"])
		if "min_covA" in k:
			params += "(GEN[ALL].DP[*] >= {}) & ".format(conf["min_covA"])
		if "min_reads_strand" in k:
			params += "((NR >= {}) & (NF >= {})) & ".format(conf["min_reads_strand"], conf["min_reads_strand"])
		if "min_reads_alt" in k:
			params += "(GEN[*].NV[*] >= {}) & ".format(conf["min_reads_alt"])
	'''elif typ == "b":
		params += "(FILTER != germline_risk) & "
		if "min_covB" in k:
			params += "(GEN[ALL].DP[*] >= {}) & ".format(conf["min_covB"])
		if "max_altB" in k:
			params += "(GEN[*].NV[*] <= {}) & ".format(conf["max_altB"])
		if "max_prop_altB" in k:
			params += "(GEN[*].NV[*] <= {}) & ".format(conf["max_prop_altB"])
	elif typ == "n":'''

	# Replace trailing spaces and ampersand with close parentheses
	params = params[:-3] + ")"
	return params	

def snpsiftFilter(conf, vcf, typ):
	# Calls bcftools to filter calls before calling isec
	if "unfiltered" in vcf:
		outfile = vcf.replace("unfiltered", "noGermline")
	log = outfile.replace("vcf", "stdout")
	if "snpsift" in conf.keys():
		cmd = ("java -jar {} filter ").format(conf["snpsift"])
	else:
		cmd = "snpsift filter "
	cmd += snpsiftOpt(conf, typ)
	cmd += " -f {} > {}".format(vcf, outfile)
	res = runProc(cmd, log)
	if res == False:
		outfile = None
	return outfile

def filterCalls(conf, vcf, typ, outdir = None):
	# Calls gatk to filter mutect calls to remove germline variants
	ext = ".unfiltered.vcf"
	if outdir:
		outfile = outdir + vcf[vcf.rfind("/")+1:vcf.find(".")] + ext
	else:
		outfile = vcf[:vcf.find(".")] + ext
	log = outfile.replace("vcf", "stdout")
	# Assemble command
	if "gatk" in conf.keys():
		cmd = ("java -jar {} FilterMutectCalls ").format(conf["gatk"])
	else:
		cmd = "gatk FilterMutectCalls "
	cmd += ("-V {} -O {}").format(vcf, outfile)
	res = runProc(cmd, log)
	if res == True and getStatus(log) == True:
		return snpsiftFilter(conf, outfile, typ)
	else:
		return None
