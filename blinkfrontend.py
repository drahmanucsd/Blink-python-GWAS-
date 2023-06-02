from haptools import data
from haptools.sim_phenotype import Haplotype
import argparse
import pandas as pd
# import statsmodels.api as sm
# import matplotlib.pyplot as plt
import os

def vcf_to_hap(vcf_path:str, hpath: str, hap_out_path: str, pheno_out_path):
    if hpath == None or hpath == '' or not os.path.isdir(hpath):
        hpath = 'haptools'
    # which variants do we want to write to the haplotype file?
    variants = {"rs149635655", "rs141306699"}

    # load the genotypes file
    # you can use either a VCF or PGEN file
    gt = data.GenotypesVCF(vcf_path)
    gt.read(variants=variants)

    # initialize an empty haplotype file
    hp = data.Haplotypes(hap_out_path, haplotype=Haplotype)
    hp.data = {}

    for variant in gt.variants:
        ID, chrom, pos, alleles = variant[["id", "chrom", "pos", "alleles"]]
        end = pos + len(alleles[1])

        # create a haplotype line in the .hap file
        # you should fill out "beta" with your own value
        hp.data[ID] = Haplotype(chrom=chrom, start=pos, end=end, id=ID, beta=0.5)

        # create variant lines for each haplotype
        hp.data[ID].variants = (data.Variant(start=pos, end=end, id=ID, allele=alleles[1]),)

    hp.write()

    #Generates the phenotype using haptools linux command
    pheno_gen_cmd = hpath + " simphenotype "+ vcf_path + " " + hap_out_path+ " -o " +pheno_out_path
    os.system(pheno_gen_cmd)

    #Reformatting the simulated phenotypes to .phen format
    df = pd.read_csv(pheno_out_path, sep='\t')
    df.insert(1, 'New Column', df.iloc[:, 0])
    df.to_csv(pheno_out_path, sep='\t', index=False)
    print(f"Simulated phenotypes saved to: {pheno_out_path}")

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
    simdata.add_argument('--hapout',type=str,help="Path to output .hap file", required=True)
    simdata.add_argument('--phenout',type=str,help="Path to output file of phenotypes", required=True)
    simdata.add_argument('--hpath', help='Haptools path (Eg: ~/.local/bin/haptools)', type=str, required=False)
    # gwas func parameters
    # takes required geno, pheno, and output file
    gwas.add_argument('--g',type=str,help="Path to the genotype file (.vcf.gz format)", required=True)
    gwas.add_argument('--p',type=str,help="Path to the test file (.phen format)", required=True)
    gwas.add_argument('--o',type=str,help="Path to output graph file (.png format)", required=True)

    args = parser.parse_args()
    
    #Function Exec calls
    if args.command == 'simdata':
        vcf_to_hap(args.i,args.hpath,args.hapout, args.phenout)
    elif args.command == 'gwas':
        abhi_cmd(args.g,args.p,args.o)
    else:
        parser.print_help()
        
main()