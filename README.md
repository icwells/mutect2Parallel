# Call Mutect2 in serial over input files
### This script is meant to meet the specific needs of my lab, but is provided in case it proves useful to others.

Copyright 2018 by Shawn Rupp

## Installation
### Download the repository:

git clone https://github.com/icwells/mutect2Serial.git

### Download dependencies:

GATK: https://software.broadinstitute.org/gatk/download/ 

Picard: https://github.com/broadinstitute/picard/releases/ 

pysam: 

	conda install pysam 

## Input files:

### Config File 
Change the name of example_config.txt (so it won't be replaced if there is another pull). 
Simply add the paths to the required files after the equals sign. 
The GATK and Picard jars can be omitted if you are loading modules on a grid system. 

### Manifest file 
The manifest file may be a space, comma, or tab seperated text file with one entry per line. 
Each entry should ahve the following format: 

	SampleID	path_to_normals	path_to_tumor_file1	path_to_tumor_file2 

### Index/Dict Files
You may include a fasta dict file and fasta and vcf indeces if they are available. 
If they are not present, the script will generate them for each file. 

## Example Usage

	python mutect2Serial.py {--jar} -i path_to_manifest -c path_to_config_file -o path_to_output_directory

Be sure that samtools is in in your PATH before running.
