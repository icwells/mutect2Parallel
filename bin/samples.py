'''This script defines classes for Samples to manage filtering of mutect2 output'''

import os
from shutil import copy
from unixpath import *
from commonUtil import *
from sample import *

class Samples():
	# Stores data for all samples in a comparison
	def __init__(self):
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
		indir = checkDir(indir)
		outdir = checkDir(outdir, True)
		self.ID = getParent(indir)
		self.Outdir = outdir + self.ID + "/"
		self.Log = self.Outdir + "mutectLog.txt"
		if not os.path.isfile(self.Log):
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
			self.appendLog(self.A)
			self.appendLog(self.B)

	def rmGermline(self):
		# Wraps calls to Sample.rmGermline
		self.A.rmGermline(self.Conf, self.Outdir)
		self.appendLog(self.A)
		self.B.rmGermline(self.Conf, self.Outdir)
		self.appendLog(self.B)

	def covB(self):
		# Calls bcf isec and uses output to generate bed files for each comparison
		run = False
		if self.A.Step == "filtering_germline" and self.A.Status == "complete":
			run = True
		elif self.A.Step == "filtering_covB" and self.A.Status != "complete":
			run = True
		if run == True:
			# Update statuses and get output file names and log
			self.updateStatuses("starting", "filtering_covB")
			self.A.Output = self.A.Unfiltered.replace(".noGermline", ".covB")
			self.B.Output = self.B.Unfiltered.replace(".noGermline", ".covB")
			# Call covB.sh: vcf1 vcf2 outputvcf2 outputvcf1 bam1 bam2 genome gatkjar
			cmd = ("bash covB.sh {} {} {} {} ").format(self.A.Private, self.B.Private, self.B.Output, self.A.Output)
			cmd += ("{} {} {} {}").format(self.A.Bam, self.B.Bam, self.Conf["ref"], self.Conf["gatk"])
			print(cmd)
			quit()
			res = runProc(cmd, log)
			if res == True:
				# Output names have already been updated
				self.updateStatuses("complete", append = True)
			else:
				self.updateStatuses("failed", append = True)

#-----------------------------------------------------------------------------

	def comparePair(self, outpath, vcfs):
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
			a = self.comparePair(aout, [self.A.Output, self.B.Unfiltered])
		else:
			a = 0
		if getTotal(self.B.Output) > 0:
			b = self.comparePair(bout, [self.B.Output, self.A.Unfiltered])
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
