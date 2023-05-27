![python](https://img.shields.io/badge/python-3.9.5-green)
![numpy](https://img.shields.io/badge/numpy-1.21.1-blue)
![pandas](https://img.shields.io/badge/pandas-1.5.3-white)
![matplotlib](https://img.shields.io/badge/matplotlib-3.4.2-red)
![project_status](https://img.shields.io/badge/project__status-work%20in%20progress-orange)


# Blink-python-GWAS
This tool, Blink, is used to perform genome-wide association studies (GWAS). GWAS is important as it helps scientists identify genes associated with diseases or traits. Our implementation will compared to a GWAS tool called plink which is well-known in the bioinformatic field.

## Installation instructions
Installation requires the `haptools` library to be installed. You can install these with pip:

`pip install haptools`

Once the required libraries are installed you can install `blink` with the following command:

`python setup.py install` (to be implemented this coming week)

Note: If you do not have root access, you can run the commands aboce with the additional options to install locally:

`pip install --user haptools` or (`pip install haptools` --> `ls ~/.local/bin/` --> `export PATH=$PATH:$HOME/.local/bin`)

`python setup.py install --user`

## Basic usage instructions


## Complete usage instructions 
We have no fully implemented the complete usage of `blink` yet and will be incorporating it this upcoming week.
* `-g FILE`, `--genotypes FILE` --> input a genotype file (.vcf.gz format or else an error will be thrown)
* `-h FILE`, `--haplotypes FILE` --> input a haplotype file (.hap format or else an error will be thrown)
* `-o FILE`, `--output FILE` --> specify an output file (for our current code a .pheno output file should be specify)

## Credits
This repository was generated by Abhishek Ganga, Daniyal Rahman, and Jung Tzen Liew with inspiration from plink, haptools, and the CSE 185 Lab 3 projects.
Please submit a pull request with any corrections or suggestions
