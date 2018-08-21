'''This script contains functions for parsing and analyzing bam/vcf files'''

import os
import gzip
import pysam
from subprocess import Popen
from shlex import split
from unixpath import getFileName

def runProc(cmd, log = None):
	# Wraps call to Popen, writes stdout/stdout err to log/devnull, returns True if no errors
	if not log:
		log = os.devnull
	with open(log, "w") as out:
		try:
			call = Popen(split(cmd), stdout = out, stderr = out)
			call.communicate()
			return True
		except:
			s = cmd.split()
			proc = s[0]
			if "-" not in s[1]:
				# Get subcommand if present
				proc += " " + s[1]
			elif proc == "java":
				# Replace call to jar with name of jar
				proc = getFileName(s[2])
			print(("\t[Warning] Could not call {}").format(proc))
			return False

def tabix(vcf, force = False):
	# tabix index and bgzips vcf files
	if os.path.isfile(vcf + ".gz") and force == False:
		gz = vcf + ".gz"
	else:
		try:
			gz = pysam.tabix_index(vcf, seq_col=0, start_col=1, end_col=1, force=True)
		except OSError:
			print(("\t[Warning] Could not index {}.").format(vcf))
			gz = None
	return gz

def bcfMerge(outfile, com):
	# Calls bftools concat on input files
	cmd = ("bcftools merge --force-samples -O v -o {} {} {}").format(outfile, com[0], com[1])
	res = runProc(cmd)
	if res == False:
		outfile = None
	return outfile	

def bcfSort(infile):
	# Call bcftools sort
	fmt = "v"
	if not infile:
		return None
	outfile = infile.replace(".vcf", ".sorted.vcf")
	if ".gz" in infile:
		fmt = "z"
	cmd = ("bcftools sort -O {} -o {} {}").format(fmt, outfile, infile)
	res = runProc(cmd)
	if res == False or os.path.isfile(outfile) == False:
		return None
	return tabix(outfile, True)

def getTotal(vcf):
	# Returns total number of content lines from vcf
	count = 0
	if os.path.isfile(vcf):
		try:
			with gzip.open(vcf, "rb") as f:
				tag = ("#").encode()
				for line in f:
					if tag not in line:
						count += 1
		except (OSError, UnicodeDecodeError):
			with open(vcf, "r") as f:
				for line in f:
					if line[0] != "#":
						count += 1
	return count

def bcfIsec(outpath, vcfs):
	# Calls bcftools to get intersecting rows and returns number of private A
	a = None
	for i in range(len(vcfs)):
		# Make sure there is an up-to-date index file
		vcfs[i] = tabix(vcfs[i], True)
	if None not in vcfs:
		cmd = ("bcftools isec {} {} -p {}").format(vcfs[0], vcfs[1], outpath)
		res = runProc(cmd)
		if res == True:
			# Number of unique variants to each sample and number of shared
			a = getTotal(outpath + "/0000.vcf")
	return a

#-----------------------------------------------------------------------------

def unifiedGenotyper(conf, samples, a, b):
	# Returns summary of unified genotyper output
	gt = samples[a].Output.replace("noGermline", "covB")
	log = gt.replace("vcf", "stdout")
	outfile = gt.replace("vcf", "tsv")
	if "gatk" in conf.keys():
		cmd = ("java -Xms512m -Xmx6G -jar {} -T UnifiedGenotyper ").format(conf["gatk"])
	else:
		cmd = "gatk -T UnifiedGenotyper "
	cmd += ("-R {} -I {} -o {} ").format(conf["ref"], samples[a].Bam, gt)
	cmd += ("--intervals {} --output_mode EMIT_ALL_SITES > {} 2>&1").format(samples[b].Bed, log)
	# Log is supplied in call to gatk, so it is not given below
	res1 = runProc(cmd)
	if res1 == True:
		sub = ('cat {} | sed "/^#/d" | ').format(gt)
		sub += "perl -lane '$F[9]=~s/^[^:]*:([^:]*).*/$1/;@reads=split(\",\",$F[9]);$reads[1]==\"\" and $reads[1]=0;if($reads[0] eq \"./.\"){$readsref=0;$readsout=0}else{$readsref=splice(@reads,0,1);$readsout=join(\",\",@reads)};print join(\"\t\",@F[0,1,3,4],$readsref,$readsout)'"
		sub += (" > {}").format(outfile)
		res2 = runProc(sub)
		if res2 == True:
			return outfile
	return None

def vcf2bed(vcf):
	# Converts target vcf to bed file, returns file name
	path = vcf[:vcf.rfind(".")]
	deletions = path + "_deletions.bed"
	snvs =  path + "_snvs.bed"
	bed = path + ".bed"
	# Get bed of deletions
	dlt = ("vcf2bed --deletions < {} > {}").format(vcf, deletions)
	res1 = runProc(dlt)
	if res1 == True:
		# Get bed of single nucleotide variants
		snv = ("vcf2bed --snvs < {} > {}").format(vcf, snvs)
		res2 = runProc(snv)
		if res2 == True:
			# Merge beds
			bdp = ("bedops --everything {} {} | awk 'BEGIN{OFS=\"\t\"}{print($1,$2,$3)}' > {}").format(deletions, snvs, bed)
			res3 = runProc(bdp)
			if res3 == True:
				return bed
	return None

#-----------------------------------------------------------------------------

def getFastaIndex(ref):
	# Creates fasta index for reference genome
	print("\tGenerating reference fasta index...")
	pysam.faidx(ref)

def samIndex(bam):
	# Call samtools to index sam/bam file
	if not os.path.isfile(bam + ".bai"):
		print(("\tGenerating sam index for {}...").format(bam))
		pysam.index(bam)

def getTumorName(bam):
	# Gets relevent header data and extracts tumor sample name from bam file
	try:
		header = pysam.view("-H", bam)
		rg = header[header.find("@RG"):header.find("@PG")]
		name = rg[rg.find("SM:"):]
		return name[name.find(":")+1:name.find("\t")]
	except pysam.utils.SamtoolsError:
		return None

def addRG(bam, sid, picard):
	# Adds readgroups to bam files
	outfile = bam[:bam.rfind(".")] + ".withRG.bam"
	if picard:
		cmd = ("java -jar {} AddOrReplaceReadGroups I={} O={} ").format(picard, bam, outfile)
	else:
		cmd = ("picard AddOrReplaceReadGroups I={} O={}").format(bam, outfile)
	cmd += (" RGLB=lib1 RGPL=illumina RGPU={}").format(sid)
	print(("\tAdding read groups to {}").format(bam))
	res = runProc(cmd)
	if res == True and outfile:
		name = getTumorName(outfile)
		return outfile, name
	else:
		return None, None

def checkRG(bam, sid, picard=None):
	# Adds read groups, creates bam index, and returns read group name
	idx = True
	name = getTumorName(bam)
	if type(name) != str:
		idx = False
		bam, name = addRG(bam, sid, picard)
		if not bam or not name:
			return None, None
		else:
			idx = True
	if idx == True:
		samIndex(bam)
		return name, bam
