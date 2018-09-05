'''This script defines classes for Samples and Sample to manage filtering of mutect2 output'''

import os
from shutil import copy
from unixpath import *
from commonUtil import *

class Samples():
	# Stores data for all samples in a comparison
	def __init__(self, parent):
		self.Ulog = ""
		self.Summary = ""
		self.Conf = {}
		self.ID = ""
		self.Outdir = ""
		self.Log = ""
		self.A = Sample()
		self.B = Sample()
		self.N = Sample()

	def appendLog(self, s):
		# Appends checkpoint status to log file
		if s.Step == "mutect" and s.Status == "starting":
			# Record infile instead of outfile
			out = s.Input
		elif "isec1" in s.Step and s.Private:
			# Record private vcf
			out = s.Private
		else:
			out = s.Output
		with open(self.Log, "a") as l:
				l.write(("{}\t{}\t{}\t{}\t{}\n").format(s.Name, s.ID, s.Step, s.Status, out))

	def setLogs(self, summary, ulog, conf):
		# Stores logs and config
		self.Summary = summary
		self.Ulog = ulog
		self.Conf = conf

	def checkSamples(self, samples):
		# Makes sure two samples have passed mutect
		proceed = True
		if "N" in samples.keys():
			self.N = samples["N"]
		else:
			printError(("Could not find normal data for {}").format(self.ID))
			proceed =  False
		if proceed == True and "A" in samples.keys():
			self.A = samples["A"]
		else:
			printError(("Could not find tumor A data for {}").format(self.ID))
			proceed =  False
		if proceed == True and "B" in samples.keys():
			self.B = samples["B"]
		else:
			printError(("Could not find tumor B data for {}").format(self.ID))
			proceed =  False
		# Ensure both have at least passed mutect
		if proceed == True:
			proceed = self.A.checkStatus(self.ID)
		if proceed == True:
			proceed = self.B.checkStatus(self.ID)				
		return proceed

	def setSamples(self, indir, outdir, done):
		# Sets samples A, B, and N; returns True if ID not in done
		ret = False
		self.ID = getParent(indir)
		self.Outdir = checkDir(checkDir(outdir, True) + self.ID)
		self.Log = outdir + "mutectLog.txt"
		if indir != outdir:
			# Copy log file to new directory
			copy(indir + "mutectLog.txt", self.Log)
		if os.path.isfile(self.Log) and self.ID not in done:
			# Proceed if sample not done and log is present
			_, s = checkOutput(self.Outdir, prnt = False)
			ret = self.checkSamples(s)
		return True

	def updateStatuses(self, status, step = None, append = False):
		# Updates A and B, appends to log if append == True
		self.A.updateStatus(status, step)
		self.B.updateStatus(status, step)
		if append == True:
			appendLog(self.Conf, self.A)
			appendLog(self.Conf, self.B)

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

	def compareVCFs(self, filt = False):
		# Compares unfilted vs. passed results for each combination of pair of samples
		outpath = conf["outpath"] + name + "/"
		if filt == False:
			aout = self.Outdir + "A_unfiltered"
			bout = self.Outdir + "B_unfiltered"
			cout = self.Outdir + "common_unfiltered.vcf"
		else:
			aout = self.Outdir + "A_filtered"
			bout = self.Outdir + "B_filtered"
			cout = self.Outdir + "common_filtered.vcf"
		# Append filtered sample name when submitting
		if getTotal(self.A.Output) > 0:
			a = comparePair(aout, [self.A.Output, self.B.Unfiltered])
		else:
			a = 0
		if getTotal(self.B.Output) > 0:
			b = comparePair(bout, [self.B.Output, self.A.Unfiltered])
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
			out.write(("{},{},{},{},{},{},{:.2%}\n").format(name, self.A.ID, self.B.ID, a, b, c, sim))
		self.A.Private = aout + "/0000.vcf"
		self.B.Private = bout + "/0000.vcf"

#-----------------------------------------------------------------------------

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
				printError(("Cannot find {} from {}").format(self.Output, sid)) 
				ret = False	
		if self.Step == "mutect":
			if self.Status != "complete":
				printError(("{} from {} has not successfully completed mutect.").format(self.ID, sid)) 
				ret = False			
		elif samples[s].Step == "filtering" and samples[s].Status != "complete":
			# Re-attempt failed filtering steps
			samples[s].Step = "mutect"
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
