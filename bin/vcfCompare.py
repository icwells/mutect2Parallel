'''This script contains fucntions for filtering and comparing vcf files'''

import os
from bamUtils import *
from commonUtil import *

def printError(msg):
	# Prints formatted error message
	print(("\n\t[Error] {}. Skipping.\n").format(msg))

#-------------------------------Comparison------------------------------------

def comparePair(outpath, vcfs):
	# Calls bcftools on pair of reads if both exist
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
	outpath = conf["outpath"] + name + "/"
	if "Unfiltered" in log:
		aout = outpath + "A_unfiltered"
		bout = outpath + "B_unfiltered"
		cout = outpath + "common_unfiltered.vcf"
	else:
		aout = outpath + "A_filtered"
		bout = outpath + "B_filtered"
		cout = outpath + "common_filtered.vcf"
	# Append filtered sample name when submitting
	if getTotal(samples["A"].Output) > 0:
		a = comparePair(aout, [samples["A"].Output, samples["B"].Unfiltered])
	else:
		a = 0
	if getTotal(samples["B"].Output) > 0:
		b = comparePair(bout, [samples["B"].Output, samples["A"].Unfiltered])
	else:
		b = 0
	if a > 0 and b > 0:
		# Merge common variants and get total and similarity
		acom = tabix(aout + "/0002.vcf")
		bcom = tabix(bout + "/0002.vcf")

		common = bcfMerge(cout, [acom, bcom])
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
		out.write(("{},{},{},{},{},{},{:.2%}\n").format(name, samples["A"].ID, samples["B"].ID, a, b, c, sim))
	return True

#-------------------------------Filtering-------------------------------------

def fmtParameters(params):
	# Adds ampersand to parameters if needed
	if len(params) > 0:
		params += " & "
	return params

def filterOpt(conf, typ):
	# Returns string with approriate parameters
	k = conf.keys()
	params = ""
	if typ == "a":
		params += "FILTER != 'germline_risk' & QUAL == '.'"
		if "min_covA" in k:
			params = fmtParameters(params)
			params += "DP[*] >= {}".format(conf["min_covA"])
		if "min_reads_strand" in k:
			params = fmtParameters(params)
			params += "F1R2[*] >= {} & F2R1[*] >= {}".format(conf["min_reads_strand"], conf["min_reads_strand"])
	elif typ == "b":
		if "min_covB" in k:
			params = fmtParameters(params)
			params += "DP[*] >= {}".format(conf["min_covB"])
		if "max_prop_altB" in k:
			params = fmtParameters(params)
			params += "AF[*] <= {}".format(conf["max_prop_altB"])
	'''elif typ == "n":'''

	return params

def bcfFilter(conf, vcf, typ):
	# Calls bcftools to filter calls before calling isec
	fmt = "v"
	if ".gz" in vcf:
		fmt = "z"
	if "unfiltered" in vcf:
		outfile = vcf.replace("unfiltered", "noGermline")
	# Get command with output format and file
	cmd = ("bcftools filter -O {} -o {} ").format(fmt, outfile)
	opt = filterOpt(conf, typ)
	if opt:
		cmd += ('-i "{}" ').format(opt)
	res = runProc(cmd + vcf)
	if res == False:
		outfile = None
	return outfile

def reheader(infile):
	# Inserts contig lines into vcf header and returns outfile name
	insert = '##INFO=<ID=P_CONTAM,Number=A,Type=Float,Description="Posterior probability an site reperesents contamination">\n'
	pos = False
	ins = True
	outfile = infile + "~"
	# Get list of ids in file and sort
	with open(outfile, "w") as out:
		with open(infile, "r") as f:
			for line in f:
				if "##INFO=" in line:
					if "P_CONTAM" in line:
						# Only insert line once
						ins = False
						break
					pos = True
				else:
					if ins == True and pos == True:
						# Insert line after last info field
						out.write(insert)
						ins = False
				out.write(line)
	# Remove infile and rename outfile
	if ins == True:
		os.remove(infile)
		os.rename(outfile, infile)
	else:
		os.remove(outfile)

def filterCalls(conf, vcf, typ, outdir = None):
	# Calls gatk to filter mutect calls to remove germline variants
	if typ == "a" and ".gz" not in vcf:
		reheader(vcf)
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
		return bcfFilter(conf, outfile, typ)
	else:
		return None
