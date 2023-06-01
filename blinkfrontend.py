# from haptools import data
# from haptools.sim_phenotype import Haplotype
import argparse
# import pandas as pd
# import statsmodels.api as sm
# import matplotlib.pyplot as plt
import os

def vcf_to_hap(vcf_path:str, hpath: str):
    if hpath == None or hpath == '' or not os.path.isdir(hpath):
        hpath = 'haptools'
    print(vcf_path,hpath)

def abhi_cmd(geno_path, pheno_path,out_path):
    print("hi")

def main():
    parser = argparse.ArgumentParser(prog='blink',description="Command-line tool to perform basic gwas")

    #Create the 2 command options
    subparser = parser.add_subparsers(dest='command')
    simdata = subparser.add_parser('simdata', help='Generates simulated phenotypes')
    gwas = subparser.add_parser('gwas', help='Generates QQ plot after perfomring gwas')

    #simdata func parameters
    #required input vcf file, and optional haptools path if default fails 
    simdata.add_argument('--i',help="Path to input vcf file", type=str, required=True)
    simdata.add_argument('--hpath', help='Haptools path (Eg: ~/.local/bin/haptools)', type=str, required=False)
    # gwas func parameters
    # takes required geno, pheno, and output file
    gwas.add_argument('--g',type=str,help="Path to the genotype file (.vcf.gz format)", required=True)
    gwas.add_argument('--p',type=str,help="Path to the test file (.phen format)", required=True)
    gwas.add_argument('--o',type=str,help="Path to output graph file (.png format)", required=True)

    args = parser.parse_args()
    
    #Function Exec calls
    if args.command == 'simdata':
        vcf_to_hap(args.i,args.hpath)
    elif args.command == 'gwas':
        abhi_cmd(args.g,args.p,args.o)
    else:
        parser.print_help()
        
main()