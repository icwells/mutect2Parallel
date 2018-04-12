'''This script contains functions for adding readgroups to bam files, as well as indexing and extracting readgroups'''

import os
import pysam
from subprocess import Popen
from shlex import split

def bgzip(vcf):
	# bgzip compresses filtered vcf files
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	H_SL-20-LK38-020-node-A1	H_SL-20-LK40-020-B3
	if os.path.isfile(vcf + ".gz"):
		gz = vcf + ".gz"
	else:	
		gz = pysam.tabix_index(vcf, seq_col=0, start_col=1, end_col=1)
	return gz

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
	if not name or len(name) < 5:
		idx = False
		bam, name = addRG(bam, sid, picard)
		if not bam or not name:
			return None, None
		else:
			idx = True
	if idx == True:
		samIndex(bam)
		return name, bam
