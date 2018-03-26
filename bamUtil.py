'''This script contains functions for adding readgroups to bam files, as well as indexing and extracting readgroups'''

import os
import pysam
from subprocess import Popen
from shlex import split

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
	# Gets relevent header data and extract tumor sample name from bam file
	header = pysam.view("-H", bam)
	rg = header[header.find("@RG"):header.find("@PG")]
	name = rg[rg.find("SM:"):]
	return name[name.find(":")+1:name.find("\t")]

def addRG(bam, sid, picard):
	# Adds readgroups to bam files
	outfile = bam[:bam.rfind(".")] + ".withRG.bam"
	if picard:
		cmd = ("java -jar {} AddOrReplaceReadGroups I={} O={} ").format(picard, bam, outfile)
	else:
		cmd = ("picard AddOrReplaceReadGroups I={} O={}").format(bam, outfile)
	cmd += (" RGLB=lib1 RGPL=illumina RGPU={}").format(sid)
	print(cmd)
	try:
		print(("\tAdding read groups to {}").format(bam))
		with open(os.devnull, "w") as dn:
			arg = Popen(split(cmd), stdout = dn, stderr = dn)
			arg.communicate()
		return outfile
	except:
		print(("\t[Error] Could not add read groups to {}\n").format(bam))
		quit()

def checkRG(bam, sid, picard=""):
	# Adds read groups, creates bam index, and returns read group name
	name = getTumorName(bam)
	if not name or len(name) < 5:
		bam = addRG(bam, sid, picard)
	samIndex(bam)
	return name, bam
