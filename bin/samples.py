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
		self.Blog = ""
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
		elif s.Step == "isec1" and s.Private:
			# Record private vcf
			out = s.Private
		else:
			out = s.Output
		with open(self.Log, "a") as l:
				l.write(("{}\t{}\t{}\t{}\t{}\n").format(s.Name, s.ID, s.Step, s.Status, out))

	def setLogs(self, summary, ulog, blog, conf):
		# Stores logs and config
		self.Summary = summary
		self.Ulog = ulog
		self.Blog = blog
		self.Conf = conf

	def __checkSamples__(self, samples):
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
			ret = self.__checkSamples__(s)
		return True

	def updateStatuses(self, status, step = None, append = False):
		# Updates A and B, appends to log if append == True
		self.A.updateStatus(status, step)
		self.B.updateStatus(status, step)
		if append == True:
			self.appendLog(self.A)
			self.appendLog(self.B)

#-----------------------------------------------------------------------------

	def rmGermline(self):
		# Wraps calls to Sample.rmGermline
		ret = False
		a = self.A.rmGermline(self.Conf, self.Outdir)
		if a == True:
			self.appendLog(self.A)
		b = self.B.rmGermline(self.Conf, self.Outdir)
		if b == True:
			self.appendLog(self.B)
		if a == True and b == True:
			ret = True
		return ret

	def __comparePair__(self, outpath, vcfs):
		# Calls bcftools on pair of reads if both exist
		for i in range(len(vcfs)):
			if vcfs[i] and os.path.isfile(vcfs[i]):
				# Make sure files are bgzipped
				vcfs[i] = tabix(vcfs[i])
			elif ".gz" not in vcfs[i] and os.path.isfile(vcfs[i] + ".gz"):
				vcfs[i] += ".gz"
			else:
				print(vcfs)
				printError(("Cannot find {}").format(vcfs[i]))
				return None
		# Call bftools on all results
		ret = bcfIsec(outpath, vcfs)
		return ret, vcfs[0], vcfs[1]

	def compareVCFs(self, step):
		# Compares unfilted vs. passed results for each combination of pair of samples
		# Get outputs and logs and check for completion
		if step == "a":
			aout = self.Outdir + "A_unfiltered"
			bout = self.Outdir + "B_unfiltered"
			cout = self.Outdir + "common_unfiltered.vcf"
			log = self.Ulog
			self.A.Private = aout + "/0000.vcf"
			self.B.Private = bout + "/0000.vcf"
			if self.B.Step != "filtering_germline" or self.B.Status != "complete":
				return False
		elif step == "b":
			aout = self.Outdir + "A_covb"
			bout = self.Outdir + "B_covb"
			cout = self.Outdir + "common_covb.vcf"
			log = self.Blog
			if self.B.Step != "filtering_forB" or self.B.Status != "complete":
				return False
		elif step == "n":
			aout = self.Outdir + "A_nab"
			bout = self.Outdir + "B_nab"
			cout = self.Outdir + "common_nab.vcf"
			log = self.Summary
			if self.B.Step != "filtering_NAB" or self.B.Status != "complete":
				return False
		# Make sure file names are updated if they are gzipped
		if getTotal(self.A.Output) > 0:
			a, self.A.Output, self.B.Unfiltered = self.__comparePair__(aout, [self.A.Output, self.B.Unfiltered])
			if step == "a":
				# Update B.Output if it was gzipped
				self.B.Output = self.B.Unfiltered
		else:
			a = 0
		if getTotal(self.B.Output) > 0:
			b, self.B.Output, self.A.Unfiltered = self.__comparePair__(bout, [self.B.Output, self.A.Unfiltered])
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
			out.write(("{},{},{},{},{},{},{:.2%}\n").format(self.ID, self.A.ID, self.B.ID, a, b, c, sim))
		return True

	def covB(self):
		# Calls covB.sh to generate bed files for each comparison
		run = False
		if self.A.Step == "isec1" and self.A.Status == "complete":
			run = True
		elif self.A.Step == "filtering_covB" and self.A.Status != "complete":
			run = True
		if run == True:
			# Update statuses and get output file names
			self.updateStatuses("starting", "filtering_covB")
			self.A.Bed = self.A.Private[:self.A.Private.rfind("/")] + "/A.private.tsv"
			self.B.Bed = self.B.Private[:self.B.Private.rfind("/")] + "/B.private.tsv"
			# Call covB.sh: vcf1 vcf2 outputvcf2 outputvcf1 bam1 bam2 genome gatkjar
			cmd = ("bash covB.sh {} {} {} {} ").format(self.A.Private, self.B.Private, self.A.Bed, self.B.Bed)
			cmd += ("{} {} {} {}").format(self.A.Bam, self.B.Bam, self.Conf["ref"], self.Conf["gatk"])
			res = runProc(cmd)
			if res == True:
				self.updateStatuses("complete", append = True)
			else:
				self.updateStatuses("failed", append = True)
		return run

	def __filterParams__(self, mode):
		# Returns parameters for heterAnalyzer
		params = ""
		if mode == "covb":
			opt = ["min_covB", "max_altB", "max_prop_altB"]
		elif mode == "nab":
			opt = ["min_covN", "max_freq_altN", "max_reads_altN"]
		for i in opt:
			if i in self.Conf.keys():
				params += ("--{} {} ").format(i, self.Conf[i])
		return params

	def filterForCov(self, mode):
		# Calls heterAnalyzer to filter for coverage in paired sample
		ret = False
		params = self.__filterParams__(mode)
		if mode == "covb":
			taga = "filtForB"
			tagb = "filtForA"
			beda = self.B.Bed
			bedb = self.A.Bed
		elif mode == "nab":
			taga = "NAB"
			tagb = "NAB"
			beda = self.N.Bed
			bedb = self.N.Bed
		a = self.A.filterForCoverage(mode, params, taga, beda)
		if a == True:
			self.appendLog(self.A)
		b = self.B.filterForCoverage(mode, params, tagb, bedb)
		if b == True:
			self.appendLog(self.B)
		if a == True or b == True:
			ret = True
		return ret

	def __unzipUnfiltered__(self):
		# Makes sure covN input is unzipped
		for i in [self.A.Unfiltered, self.B.Unfiltered]:
			if getExt(i) == "gz":
				res = runProc((("gzip -d {}").format(i)))
		self.A.Unfiltered = self.A.Unfiltered.replace(".gz", "")
		self.B.Unfiltered = self.B.Unfiltered.replace(".gz", "")

	def covN(self):
		# Calls covN.sh to extract coverage from normal bam file
		run = False
		if self.A.Step == "isec2" and self.A.Status == "complete":
			run = True
		elif self.N.Step == "filtering_covN" and self.N.Status != "complete":
			run = True
		if run == True:
			self.__unzipUnfiltered__()
			self.N.Bed = self.Outdir + "normalVariants.tsv"
			# Assign bed as outfile so it is recorded in log
			self.N.updateStatus("starting", "filtering_covN", self.N.Bed)
			cmd = ("bash covN.sh {} {} {} {}").format(self.A.Unfiltered, self.B.Unfiltered, self.N.Bed, self.N.Bam)
			cmd += (" {} {}").format(self.Conf["ref"], self.Conf["gatk"])
			print(cmd)
			res = runProc(cmd)
			if res == True:
				self.N.updateStatus("complete")
			else:
				self.N.updateStatus("failed")
			self.appendLog(self.N)
		return run
