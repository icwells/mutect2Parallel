'''This script contains functions for parsing and analyzing bam/vcf files'''

import os
import pysam
from subprocess import Popen
from shlex import split

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

def bcfMerge(path, com):
	# Calls bftools concat on input files
	outfile = path + "common.vcf"
	cmd = ("bcftools merge --force-samples -O v -o {} {} {}").format(outfile, com[0], com[1])
	with open(os.devnull, "w") as dn:
		try:
			bc = Popen(split(cmd), stdout = dn, stderr = dn)
			bc.communicate()
		except:
			print(("\t[Error] calling bcftools merge on samples in {}").format(path))
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
	with open(os.devnull, "w") as dn:
		try:
			bs = Popen(split(cmd), stdout = dn, stderr = dn)
			bs.communicate()
		except:
			print(("\t[Error] calling bcftools sort on {}").format(infile))
			return None
	if not os.path.isfile(outfile):
		return None
	return tabix(outfile, True)

def getTotal(vcf):
	# Returns total number of content lines from vcf
	count = 0
	if os.path.isfile(vcf):
		with open(vcf, "r") as f:
			for line in f:
				if line[0] != "#":
					count += 1
	return count

def bcfIsec(outpath, vcfs):
	# Calls bcftools to get intersecting rows and returns number of private A
	for i in range(len(vcfs)):
		# MAke sure there is an up-to-date index file
		vcfs[i] = tabix(vcfs[i], True)
	if None in vcfs:
		return None
	cmd = ("bcftools isec {} {} -p {}").format(vcfs[0], vcfs[1], outpath)
	try:
		bcf = Popen(split(cmd))
		bcf.communicate()
		summarize = True
	except:
		print(("\t[Error] Could not call bcftools isec with {}").format(cmd))
		return None
	# Number of unique variants to each sample and number of shared
	a = getTotal(outpath + "/0000.vcf")
	return a

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
	try:
		print(("\tAdding read groups to {}").format(bam))
		with open(os.devnull, "w") as dn:
			arg = Popen(split(cmd), stdout = dn, stderr = dn)
			arg.communicate()
	except:
		print(("\t[Error] Could not add read groups to {}\n").format(bam))
		return None, None
	if outfile:
		name = getTumorName(outfile)
		return outfile, name

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
