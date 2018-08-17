# mutect2Parallel v0.3 will create batch scripts to call Mutect2 in parallel over input files

Copyright 2018 by Shawn Rupp

## Installation
### Download the repository:

git clone https://github.com/icwells/mutect2Parallel.git

### Download dependencies:
#### Make sure bcftools and bedops are in your PATH

GATK: https://software.broadinstitute.org/gatk/download/ 
Picard: https://github.com/broadinstitute/picard/releases/ 
snpeff: http://snpeff.sourceforge.net/download.html  

bedops: https://bedops.readthedocs.io/en/latest/index.html 
bcftools: http://www.htslib.org/download/  

pysam: 

	conda install pysam  

unixpath:  

	git clone https://github.com/icwells/unixpath.git    
	cd unixpath/  
	python setup.py install 

## Input files:

### Config File 
Change the name of example_config.txt (so it won't be replaced if there is another pull). 
Simply add the paths to the required files after the equals sign. 
The GATK and Picard jars can be omitted if you are loading modules on a grid system. 
Be sure to include a template batch file after the variable definitions. The appropriate command 
will be added to the template and the job name will have the respective sample ID appended to it. It is
recommended that mutect2 has at least 16Gb of RAM available. Since each batch script will call two parallel 
instances, it is best to specify 32Gb of RAM for each script. 

	reference_genome		Path to reference genome (required). 
	bed_annotation			Path to BED annotation (greatly speeds up runtime) 
	output_directory		Path to parent output directory (each sample will have it's own folder) 
	GATK_jar	 			Path to GATK jar (omit if using a module). 
	Picard_jar	 			Path to Picard jar (omit if using a module). 
	normal_panel	 		Path to panel of normals VCF (Can use mutect2Parallel/getPON to generate)
	germline_resource		Path to germline recource file. 
	allele_frequency		Decimal frequency of alleles in germline recource (i.e. 1/# of individuals; required if using germline resource). 
	contaminant_estimate	Path to VCF used to estimate contmaination in samples. 

The followng lines are to include any additional Mutect or FilterMutect options. Simply enter the flag and option as you would for 
gatk (i.e. --enable_clustered_read_position_filter --annotation {annotator}). Be aware that these option will not be checked for errors by 
the python script, so they may lead to errors when calling gatk if they are not correct. 

	mutect_options			Additional options for mutect2
	filter_mutect_options	Additional options for FilterMutectCalls

### Manifest file 
The manifest file may be a space, comma, or tab seperated text file with one entry per line. 
Each entry should have the following format: 

	SampleID	path/to/normals	path/to/tumor/file1	path/to/tumor/file2 

The same format can be used to generate a new panel of normals, but only the normal file will be used. 

### Index/Dict Files
You may include a fasta dict file and fasta, bam, and vcf indeces if they are available. 
If they are not present, the script will generate them for each file. 

## Example Usage
mutect2Parallel will create one batch script per sample in the manifest file using the template provided in the config file. 
It will also check for indexes for any reference files (i.e. a reference fasta index) and generate them if necessary. 
The resulting batch scripts will run each tumor-normal combination in parallel for each sample. 

	python mutect2Parallel.py {--submit/bamout/newPON} -i path/to/manifest -c path/to/config/file -o path/to/output/directory

	-h, --help		show this help message and exit
	--submit		Submit batch files to SLURM/Torque grid for execution.
	--bamout		Indicates that mutect should also generate bam output files (extends mutect runtime).
	--newPON		Creates batch scripts for running mutect in tumor-only mode on normals 
						and creating a panel of normals (instead of running both tumor-normal comparisons)
	-i I			Path to space/tab/comma seperated text file of input files (format: ID Normal A B)
	-c C			Path to config file containing reference genome, java jars (if using), and mutect options.
	-o O			Path to batch script output directory (leave blank for current directory).

After all of the batch scripts have finished running filterVCFs.py can be used to filter the mutect output and compare the resulting vcfs 
using bcftools isec. Each filtered file will be compared to the unfiltered vcf of the other sample (i.e. filtered A vs unfiltered B 
and vice versa) first using default parameters and then using the "-f .,PASS" options. 

	python filterVCFs.py {--summarize} -c path/to/config/file 

	-h, --help		show this help message and exit
	-c C			Path to config file containing reference genome, java jars (if using), and mutect options 
						(required; input files are read from sub-directories in output_directory and output will be written to same sub-directory).  
	-o				Option output dirctory to write filtered vcf files to. It will have the same structure as the mutect output, but in a 
						seperate direcotry to avoid overwriting other filtering output.  
	-t				Number of threads.  

## Other Scripts
runPair and getPON commands are formatted in batch scripts by mutect2Parallel, so it may not be necessary to directly call either. 

### runPair.py
Used to call mutect2 in parallel for each each tumor-normal comparison for one sample. This script is called by mutect2Parallel.py by default. 

	-h, --help		show this help message and exit
	--bamout		Indicates that mutect should also generate bam output files.
	-s S			Sample name (required).
	-x X			Path to first tumor bam (required).
	-y Y			Path to second tumor bam (required).
	-c C			Path to normal/control bam (required).
	-r R			Path to reference genome (required).
	-o O			Path to output directory (required).
	--bed BED		Path to bed annotation.
	--gatk GATK		Path to gatk jar (if using).
	--picard PICARD		Path to picard jar (if using).
	-p P			Path to panel of normals.
	-g G			Path to germline resource.
	--af AF			Estimated allele frequency (required if using a germline resource).
	--mo MO			Additional mutect options in quotes


### compareVariants.py  
This script will compare variants from different filtering pipelines.

	-h, --help show this help message and exit
	-c C		Copy target platypus data to this directory.
	-v V		Path to uncompressed vcf header (Copies contig information to
					platypus vcf headers).
	-i I		Path to input manifest (mutect input for manifest generation or
					generated manifest for comparison).
	-m M		Path to mutect2parallel parent output directory.
	-p P		Path to platypus-based parent output directory.
	-o O		Path to output manifest if using -m and -p. Path to output
					directory if using -i.

### getPON.py
Can be used to genrate a new panel of normals. This script will be called by mutect2Parallel.py if the --newPON flag is given. 

	-h, --help		show this help message and exit
	--pon			Generate new panel of normals from log file (requires -l
						(output from tumor only mode) and -o (output PON file) flags only).
	-s S			Sample name (required).
	-l L			Path to log file (required; output files are recorded here).
	-c C			Path to normal/control bam (required).
	-r R			Path to reference genome (required).
	-o O			Path to output directory (required).
	--bed BED		Path to bed annotation.
	--gatk GATK		Path to gatk jar (if using).
	--picard PICARD		Path to picard jar (if using).
	-g G			Path to germline resource.
	--af AF			Estimated allele frequency (required if using a germline resource).
	-e E			Path to contmination estimate vcf.


### getActiveRegion.py 
Can be used to subset a bed annotation to examine specific regions. 

	-h, --help		show this help message and exit
	-c C			Chromosome(s) to subset (seperate with commas if there is more than one).
	-i I			Path to input file.
	-o O			Path to output file.

