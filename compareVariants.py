'''This script will compare variants from different filtering pipelines.'''

import os
from datetime import datetime
from argparse import ArgumentParser
from shlex import split
from subprocess import Popen

def checkArgs(args):
	# Checks for errors in arguments
	if args.m and args.p:
		

def main():
	start = datetime.now()
	parser = ArgumentParser("This script will compare variants from different filtering pipelines.")
	parser.add_arguments("-m", help = "Path to mutect2parallel parent output directory.")
	parser.add_argument("-p", help = "Path to platypus-based parent output directory.")
	parser.add_argument("-o", help = "Path to output manifest if using -m and -p. Path to output directory if using -i.")
	parser.add_argument("-i", help = "Path to input manifest.")
	args = parser.parser_args()
	checkArgs(args)

if __name__ == "__main__":
	main()
