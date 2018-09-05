'''This script defines classes for Sample to manage filtering of mutect2 output'''

import os
from unixpath import *
import commonUtil

class Sample():
	# Stores data for managing sample progress
	def __init__(self):
		self.Name = ""
		self.ID = ""
		self.Step = ""
		self.Status = ""
		self.Output = ""
		self.Private = ""
		self.Bed = ""
		self.Bam = ""
		self.Input = ""
		self.Unfiltered = ""

	def __str__(self):
		# Return formatted string
		ret = "Name: {}\n".format(self.Name)
		ret += "ID: {}\n".format(self.ID)
		ret += "Step: {}\n".format(self.Step)
		ret += "Status: {}\n".format(self.Status)
		ret += "Output VCF: {}\n".format(self.Output)
		ret += "Private Variants: {}\n".format(self.Private)
		ret += "Bed File: {}\n".format(self.Bed)
		ret += "Source Bam File: {}\n".format(self.Bam)
		ret += "Input VCF: {}\n".format(self.Input)
		ret += "Unfiltered VCF: {}\n".format(self.Unfiltered)
		return ret

	def update(self, sample, name, step, status, outfile):
		# Sorts and updates entry with additional status update
		save = False
		if not self.Name:
			self.Name = sample
		if not self.ID:
			self.ID = name
		if step == "comparison":
			save = True
		elif step == "filtering_covB" and self.Step != "comparison":
			save = True
		elif step == "filtering_germline" and self.Step == "mutect":
			unf = outfile[:outfile.find(".")] + ".noGermline.vcf"
			if os.path.isfile(unf):
				self.Unfiltered = unf
			elif os.path.isfile(unf + ".gz"):
				self.Unfiltered = unf + ".gz"
			save = True
		elif step == "mutect" and not self.Step or self.Step == "mutect":
			save = True
		elif step == "normal" and not self.Step:
			save = True
		if save == True:
			self.Step = step
			self.Output = outfile
			self.Status = status
		if getExt(outfile) == "bam":
			self.Bam = outfile

	def checkStatus(self, sid):
		# Makes sure filtering can proceed from mutect output
		ret = True
		if not os.path.isfile(self.Output):
			if os.path.isfile(self.Output + ".gz"):
				self.Output = self.Output + ".gz"
			else:
				commonUtil.printError(("Cannot find {} from {}").format(self.Output, sid)) 
				ret = False	
		if self.Step == "mutect":
			if self.Status != "complete":
				commonUtil.printError(("{} from {} has not successfully completed mutect.").format(self.ID, sid)) 
				ret = False			
		elif self.Step == "filtering" and self.Status != "complete":
			# Re-attempt failed filtering steps
			self.Step = "mutect"
		return ret

	def updateStatus(self, status, step = None, outfile = None, unfilt = False):
		# Updates current status of sample
		self.status = status
		if step:
			self.Step = step
		if outfile:
			self.Output = outfile
			if getExt(outfile) == "bam":
				self.Bam = outfile
		if unfilt == True:
			self.Unfiltered = outfile

#-------------------------------Filtering-------------------------------------

	def fmtParameters(self, params):
		# Adds ampersand to parameters if needed
		if len(params) > 0:
			params += " & "
		return params

	def filterOpt(self, conf):
		# Returns string with approriate parameters
		k = conf.keys()
		params = ""
		if self.Step == "filtering_germline":
			params += "FILTER != 'germline_risk' & QUAL == '.'"
			if "min_covA" in k:
				params = self.fmtParameters(params)
				params += "DP[*] >= {}".format(conf["min_covA"])
			if "min_reads_strand" in k:
				params = self.fmtParameters(params)
				params += "F1R2[*] >= {} & F2R1[*] >= {}".format(conf["min_reads_strand"], conf["min_reads_strand"])
		elif self.Step == "filtering_covB":
			if "min_covB" in k:
				params = self.fmtParameters(params)
				params += "DP[*] >= {}".format(conf["min_covB"])
			if "max_prop_altB" in k:
				params = self.fmtParameters(params)
				params += "AF[*] <= {}".format(conf["max_prop_altB"])
		'''elif self.Step == "filtering_NAB":'''

		return params

	def bcfFilter(self, conf):
		# Calls bcftools to filter calls before calling isec
		fmt = "v"
		if ".gz" in self.Output:
			fmt = "z"
		if self.Step == "filtering_germline":
			outfile = self.Output.replace("unfiltered", "noGermline")
		# Get command with output format and file
		cmd = ("bcftools filter -O {} -o {} ").format(fmt, outfile)
		opt = self.filterOpt(conf)
		if opt:
			cmd += ('-i "{}" ').format(opt)
		res = commonUtil.runProc(cmd + self.Output)
		# Record results
		if res == True and os.path.isfile(outfile):
			self.updateStatus("complete", outfile = outfile, unfilt = True)
		else:
			self.updateStatus("failed")

	def reheader(self):
		# Inserts contig lines into vcf header and returns outfile name
		insert = '##INFO=<ID=P_CONTAM,Number=A,Type=Float,Description="Posterior probability an site reperesents contamination">\n'
		pos = False
		ins = True
		outfile = self.Output + "~"
		# Get list of ids in file and sort
		with open(outfile, "w") as out:
			with open(self.Output, "r") as f:
				for line in f:
					if "##INFO=" in line:
						if "P_CONTAM" in line:
							# Only insert line once
							ins = False
						pos = True
					else:
						if ins == True and pos == True:
							# Insert line after last info field
							out.write(insert)
							ins = False
					out.write(line)
		# Remove infile and rename outfile
		os.remove(self.Output)
		os.rename(outfile, self.Output)

	def filterCalls(self, conf, outdir):
		# Calls gatk to filter mutect calls to remove germline variants
		if ".gz" not in self.Output:
			self.reheader()
		ext = ".unfiltered.vcf"
		outfile = outdir + self.Output[self.Output.rfind("/")+1:self.Output.find(".")] + ext
		log = outfile.replace("vcf", "stdout")
		# Assemble command
		cmd = ("java -jar {} FilterMutectCalls ").format(conf["gatk"])
		cmd += ("-V {} -O {}").format(self.Output, outfile)
		res = commonUtil.runProc(cmd, log)
		if res == True and commonUtil.getStatus(log) == True and commonUtil.getTotal(outfile) > 0:
			self.Output = outfile
			self.bcfFilter(conf)
		else:
			self.updateStatus("failed")

	def rmGermline(self, conf, outdir):
		# Calls filterMutectCalls and bcftools to remove germline risks
		run = False
		if self.Step == "mutect" and self.Status == "complete":
			run = True
		elif self.Step == "filtering_germline" and self.Status == "failed":
			run == True
		if run == True:
			self.updateStatus("starting", "filtering_germline")
			unfiltered = self.filterCalls(conf, outdir)
		return run
