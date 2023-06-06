![python](https://img.shields.io/badge/python-3.9.5-green)
![matplotlib](https://img.shields.io/badge/matplotlib-3.4.2-red)
![numpy](https://img.shields.io/badge/numpy-1.21.1-blue)
![pandas](https://img.shields.io/badge/pandas-1.5.3-white)
![project_status](https://img.shields.io/badge/project__status-work%20in%20progress-orange)

# Blink-python-GWAS
This tool, Blink, is used to perform genome-wide association studies (GWAS). GWAS is important as it helps scientists identify genes associated with diseases or traits. Our implementation will compared to a GWAS tool called plink which is well-known in the bioinformatic field. 

## Table of Contents
* Installation Instructions
* Example Testing
* Usage Instructions
* Credits

 **_NOTE:_**  All citation and resources are located in the Notes.md file

## Installation Instructions
1. Run this on your datahub server! (All the libraries are installed on there already other than haptools)
Installation requires the [`haptools`](https://haptools.readthedocs.io/en/stable/project_info/installation.html) library to be installed. You can install these with pip:

>`pip install haptools`

Note: If you do not have root access, you can run the commands aboce with the additional options to install locally:

>`pip install --user haptools` or (`pip install haptools` --> `ls ~/.local/bin/` --> `export PATH=$PATH:$HOME/.local/bin`)

2. Then git clone our repository using this command (We recommend making a new directory and then cloning it in that directory):

>`git clone https://github.com/drahmanucsd/Blink-python-GWAS-.git`

3. Once you are in the directory you created and the required libraries are installed you can install `blink` with the following command:

>`cd Blink-python-GWAS-`
`python setup.py install`

or
>`python setup.py install --user` (if you don't have root permission and run into an error) 

**_DO THIS RIGHT AFTER AS WELL_** `export PATH=$PATH:$HOME/.local/bin`

If `blink` was installed correctly, try `blink -h` to see if instructions to run `blink` are provided.

## Testing Data
### Simulating phenotypes to run GWAS on

* Download the vcf file we used using 

>`wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz`
* It will take a while because it is an over 2gbs
* Feel free to change the variants in the blink.py file but we have defaulted it to rs149635655 and rs141306699 to work with that specific 1000Genomes dataset

Run the command below to simulate the phenotypes:

>`blink simdata --i GENOTYPES.vcf.gz --hapout PATH_TO_output.hap --phenout NAME_OF_FILE.phen --hpath ~/.local/bin/haptools`

Ex: `blink simdata --i ~/teams/31/ALL.chr21.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz --hapout ~/teams/31/output.hap --phenout simulated.phen --hpath ~/.local/bin/haptools`

If you would like more information run this command: `python blink.py simdata -h`

### Running GWAS

* Make sure you have a .vcf.gz file and .phen file. Most of these are generated or retrieved from the step above

Run the command below to run GWAS on the dataset:

For this part of the test we will be using lab3_GWAS's .vcf.gz and .phen file due to the size being smaller

>`blink gwas --g GENOTYPES.vcf.gz --p NAME_OF_FILE.phen --o PATH_OF_OUTPUT`

Run this example command below to run it on the lab 3 GWAS data!

Ex: `blink gwas --g ~/public/lab3/lab3_gwas.vcf.gz --p ~/public/lab3/lab3_gwas.phen`

Access the test.png file which will be located in the directory you call the command on and will contain the Manhattan and QQ plot!

If you would like more information run this command: 

>`python blink.py gwas -h`

Have fun :)

## Basic usage instructions
Using the python/jupyter notebook files following the instructions in the cell blocks/comments to run the code. You will most likely need to fill in your corresponding files/data to run the code. FOLLOW THE INSTRUCTIONS PROVIDED!
* For the simulate_pheno jupyter notebook, you will need to use your own variants to run the tool for converting .vcf files to .hap files

The basic usage of `blink` will be: 

`blink [other options]` (will be specified soon once we integrate everything)

## Complete usage instructions 
The following are the currently working options/parameters available for the phenotype simulation and linear regression GWAS

### Options to Simulate Phenotypes (Ran using the argument `simdata` followed by the options below)
* `--i` --> input the path of genotype file (.vcf.gz format or else an error will be thrown)
* `--hapout` --> input the path of haplotype file (.hap format or else an error will be thrown)
* `--phenout` --> specify an output file for the phenotypes (for our current code a .phen output file should be specify)
* `--hpath` --> this is just the location of where haptools is installed ~/.local/bin/haptools

### Options to run GWAS (Ran using the argument `gwas` followed by the options below)
* `--geno` --> input the path of genotype file (.vcf format or else an error will be thrown)
* `--pheno` --> input the path of phenotype file (.phen format or else an error will be thrown)
* `--out` --> output file name of association pvals and png graphs (no file format because ouptuts multiple files)
* `--maf` --> minor allelic frequency threshold for filtering low data SNPs (Optional cmd, default=0.01)

## Credits
This repository was generated by Abhishek Ganga, Daniyal Rahman, and Jung Tzen Liew with inspiration from plink, haptools, and the CSE 185 Lab 3 projects. Please submit a pull request with any corrections or suggestions! We accept all criticism.
