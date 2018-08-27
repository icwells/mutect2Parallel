'''This script will check the output logs from mutect2parallel to identify 
any samples which were not successful'''

from argparse import ArgumentParser
from glob import glob
from unixpath import checkDir

def identifyFails(outdir):
	# Searches mutect logs to identify samples which did not complete successfully
	paths = glob(outdir + "*/mutectLog.txt")
	for i in paths:
		# Isolate sample names
		with open(i, "r") as f:
			lines = f.readlines()
		for j in lines[-2:]:
			if "complete" not in j:
				print(j.strip())

def main():
	parser = ArgumentParser("This script will check the output logs from \
mutect2parallel to identify any samples which were not successful.")
	parser.add_argument("outdir", help = "Path to output directory of mutect2 parallel.")
	args = parser.parse_args()
	outdir = checkDir(args.outdir)
	identifyFails(outdir)

if __name__ == "__main__":
	main()
