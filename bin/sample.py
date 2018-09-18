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

	def __rankStatus__(self, step, status, outfile):
		# Determines if new data should overwrite existing data
		save = False
		if not self.Step:
			# Save if fields are empty
			save = True
		elif step == self.Step and self.Status != "complete":
			save = True
		elif step == "isec3":
			save = True
		elif step == "filtering_NAB" and self.Step != "isec3":
			save = True
		elif step == "filtering_covN" and self.Step == "isec2":
			save = True
		elif step == "isec2" and self.Step == "filtering_forB":
			save = True
		elif step == "filtering_forB" and self.Step == "filtering_covB":
			save = True
		elif step == "filtering_covB" and self.Step == "isec1":
			save = True
		elif step == "isec1" and self.Step == "filtering_germline":
			self.Private = outfile
			save = True
		elif step == "filtering_germline" and self.Step == "mutect":
			# Save unfiltered file if it exists
			unf = outfile[:outfile.find(".")] + ".noGermline.vcf"
			if os.path.isfile(unf):
				self.Unfiltered = unf
			elif os.path.isfile(unf + ".gz"):
				self.Unfiltered = unf + ".gz"
			save = True
		elif step == "mutect" and self.Step == "normal":
			save = True
		elif step == "normal":
			save = True
		return save

	def update(self, sample, name, step, status, outfile):
		# Sorts and updates entry with additional status update
		if not self.Name:
			self.Name = sample
		if not self.ID:
			self.ID = name
		save = self.__rankStatus__(step, status, outfile)
		if save == True:
			self.Step = step
			self.Status = status
			if step != "isec1":
				self.Output = outfile
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

	def reset(self):
		# Resets status to begin filtering
		self.Step = "mutect"
		self.Status = "complete"
		self.Unfiltered = ""

#-------------------------------Filtering-------------------------------------

	def __fmtParameters__(self, params):
		# Adds ampersand to parameters if needed
		if len(params) > 0:
			params += " & "
		return params

	def __filterOpt__(self, conf):
		# Returns string with approriate parameters
		k = conf.keys()
		params = "FILTER != 'germline_risk' & QUAL == '.'"
		if "min_covA" in k:
			params = self.__fmtParameters__(params)
			params += "DP[*] >= {}".format(conf["min_covA"])
		if "min_reads_strand" in k:
			params = self.__fmtParameters__(params)
			params += "F1R2[*] >= {} & F2R1[*] >= {}".format(conf["min_reads_strand"], conf["min_reads_strand"])
		return params

	def bcfFilter(self, conf):
		# Calls bcftools to filter calls before calling isec
		fmt = "v"
		if ".gz" in self.Output:
			fmt = "z"
		outfile = self.Output.replace("unfiltered", "noGermline")
		# Get command with output format and file
		cmd = ("bcftools filter -O {} -o {} ").format(fmt, outfile)
		opt = self.__filterOpt__(conf)
		if opt:
			cmd += ('-i "{}" ').format(opt)
		res = commonUtil.runProc(cmd + self.Output)
		# Record results
		if res == True and os.path.isfile(outfile):
			self.updateStatus("complete", outfile = outfile, unfilt = True)
		else:
			self.updateStatus("failed")

	def __reheader__(self):
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
			self.__reheader__()
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
			self.filterCalls(conf, outdir)
		return run

	def filterForCoverage(self, mode, params, tag, bed):
		# Filters for coverage using given mode and parameters
		run = False
		if mode == "covb":
			step1 = "filtering_covB"
			step2 = "filtering_forB"
		elif mode == "nab":
			step1 = "isec2"
			step2 = "filtering_NAB"			
		if self.Step == step1 and self.Status == "complete":
			run = True
		elif self.Step == step2 and self.Status != "complete":
			run = True
		if run == True:
			# Update statuses and get output file names
			infile = self.Output
			self.Output = ("{}.{}.vcf").format(self.Output[:self.Output.find(".")], tag)
			self.updateStatus("starting", step2, self.Output)
			cmd = ("./heterAnalyzer {} {}").format(mode, params)
			cmd += ("-v {} -i {} -o {}").format(infile, bed, self.Output)
			res = commonUtil.runProc(cmd)
			if res == True:
				self.updateStatus("complete")
			else:
				self.updateStatuses("failed")
		return run
