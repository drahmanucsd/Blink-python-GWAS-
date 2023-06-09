from haptools import data
from haptools.sim_phenotype import Haplotype
import numpy as np
import seaborn as sns
import statsmodels.api as sm
from qqman import qqman
import argparse
import pandas as pd
import matplotlib
# import statsmodels.api as sm
from matplotlib import pyplot as plt
import os
import warnings
warnings.filterwarnings("ignore")

## Data Simulation Functions

def vcf_to_hap(vcf_path:str, hpath: str, hap_out_path: str, pheno_out_path):
    """
    Function: Generates phenotype files given vcf input
    Parameters:
      Inputs:
        vcf_path: input vcf file use to gen phenotype file
        hpath: path containing haptools program if saved in different PATH
      Outputs:
        hap_out_path: output file contains haploids to run haptools 'simphenotypes'
        pheno_out_path: output file for the simulated phenotypes

    THE CODE BELOW IS CREDITED TO:
        Massarat, A. R., Lamkin, M., Reeve, C., Williams, A. L., D’Antonio, M., & Gymrek,
        M. Haptools: a toolkit for admixture and haplotype analysis [Computer software].
        https://github.com/CAST-genomics/haptools who created haptools 
        (Was granted permission to use from Professor Gymrek)
    """
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

    # Remove the first line
    with open(pheno_out_path, 'r') as file:
        lines = file.readlines() 
    lines = lines[1:]
    # Open the file in write mode and write the modified lines
    with open(pheno_out_path, 'w') as file:
        file.writelines(lines)

    print(f"Simulated phenotypes saved to: {pheno_out_path}")


##GWAS Linear Regression Functions

def process_SNPS(genotype_df, phenotype_df, row_number,maf):
    """
    Funtion: Process and generate a beta and pval for each SNP/row of geno data
    Parameters:
        Inputs:
            genotype_df: pandas df of genomic data
            phenotype_df: pandas df of phenotype data
            row_number: global row number of SNP in genotype data file
            maf: minor allelic frequency used for filtering

    """
    #CREATE THE ARRAY OF GENOTYPE DATA
    column_names = genotype_df.columns.tolist() #this is a list of all the column names/individuals
    values = genotype_df.iloc[row_number].values.tolist()
    #print(values)
    ressive_count = 0
    dominant_count = 0
    for val in values:
        ressive_count += val.count('1')
        dominant_count += val.count('0')
    if ressive_count/(ressive_count+dominant_count)<maf:
        return None, None
        
    sum_values = [int(string.split('|')[0]) + int(string.split('|')[1][0]) for string in values] #this is a list of the summed values per row

    gts = np.array(sum_values)
    
    if np.var(gts) == 0:
        # print(f"Zero variance for row {row_number}")
        return None, None
    
    gts = (gts-np.mean(gts))/np.sqrt(np.var(gts))
    pts = phenotype_df.get("Val").to_numpy()
    
    #getting beta and pval
    X = sm.add_constant(gts)
    model = sm.OLS(pts, X)
    results = model.fit()
    beta = results.params[1]
    pval = results.pvalues[1]
    
    return beta, pval

def plot_qq(pval_list, beta_list, genotype_data,outfile):
    """
    Funtion: Generates the QQ plot and Manhattan plot side-by-side
    Parameters:
        Inputs:
            pval_list: list of the pvalues generated from 'process_SNPS'
            beta_list: list of beta values generated from 'process_SNPS'
            genotype_data: the genomic data, 0: homo-ref, 1: hetro, 2: homo-alt
        Outputs:
            outfile: output dir of both graphs in png format
            
    """
    zero_list = [0] * len(pval_list) #NEED TO CHANGE THIS
    # zero_list = [x/0.15 for x in beta_list]
    #first make the dataframe for a qq plot
    qq_df = pd.DataFrame({'CHR': genotype_data.get("#CHROM"), 'SNP': genotype_data.get("ID"), 
                      "BP" : genotype_data.get("POS"), "A1" : genotype_data.get("REF"), "TEST" : genotype_data.get("FILTER"), 
                      "NMISS" : genotype_data.get("QUAL"), "BETA" : beta_list, "STAT" : zero_list, "P" : pval_list})
    
    fig, (ax0, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
    fig.set_size_inches((15, 5))
    qqman.manhattan(qq_df, ax=ax0)
    qqman.qqplot(qq_df, ax=ax1)

    # Save assoc file and graphs with name outfile
    csv_path = outfile + '.csv'
    fig_path = outfile + '.png'
    qq_df.to_csv(csv_path)
    plt.savefig(fig_path, bbox_inches='tight') 


def find_skip_lines(file_path):
    """
    Function: Searches Genome file and outputs the terminating line of the header
    Parameter file_path: the user inputed genome file
    """
    with open(file_path, 'r') as file:
        for line_number, line in enumerate(file, 1):
            if '#CHROM' in line:
                return line_number-1
    return 257  # Default value most likely to be

def readData(genotypeData, phenotypeData,outfile,maf=0.01):
    """
    Funtion: Reads in the data from the genotype and phenotype files and converts them for processing
    Parameters:
        Inputs:
            genotypeData: path to genome dataset in uncompressed vcf format
            phenotypeData: path to phenotype dataset, in phen format
            maf: used for filtering any SNPs below minor-allelic-freq threshold
        Outputs:
            outfile: file name of association pvals and png graph generation from GWAS
    
    """
    pval_list = [] #hold pvals
    beta_list = [] #holds beta vals
    
    lines_to_skip = find_skip_lines(genotypeData)
    # Reads in genotype data
    genotypes = pd.read_csv(genotypeData, skiprows=lines_to_skip, sep='\t')
    # Reads in phenotype data
    phenotypes = pd.read_csv(phenotypeData, sep='\t', header=None, names=['ID', 'ID2', 'Val']).drop(columns=["ID2"])
    # get_snps = genotypes.head(len(genotypes))
    
    # This runs our code on a small chunk of data
    get_snps = genotypes.head(100)
    get_snps_reformat = get_snps.drop(columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
    
    # Reorganize the phenotypes df such that it lines up with the order of the samples that occur in the genotypes file
    for name_id in get_snps_reformat.columns.values.tolist():
        temp = phenotypes.loc[phenotypes['ID'] == name_id]
        if temp.empty:
            val = 0
        else:
            val = temp.iloc[0]['Val']
        phenotypes_reformat = phenotypes_reformat.append({'ID':name_id, 'Val':val},ignore_index=True);
    
    #Iterates over every row in genotypes to generate a pval and beta list for plotting
    for i in range(get_snps_reformat.shape[0]):
        out_beta, out_pval = process_SNPS(get_snps_reformat, phenotypes_reformat, i,maf) #CALL process_SNPS
        pval_list.append(out_pval)
        beta_list.append(out_beta)

    #call plot_data using pval_list and beta_list
    plot_qq(pval_list, beta_list, get_snps,outfile)
    # plt.savefig(outfile, bbox_inches='tight') #let user choose the output file name??

def main():
    """
    THE CODE BELOW WAS INSPIRED BY:
    ARGPARSE tutorial: -- Adapted the framwork for our command line options
    https://realpython.com/command-line-interfaces-python-argparse/
    """
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
    gwas.add_argument('--geno',type=str,help="Path to the genotype file (.vcf.gz format)", required=True)
    gwas.add_argument('--pheno',type=str,help="Path to the test file (.phen format)", required=True)
    gwas.add_argument('--out',type=str,help="Path to output graph file (.png format)", required=True)
    gwas.add_argument('--maf',type=str,help="Removes SNPS with too much genotype data missing based on minor allelic frequency (defualt 0.01)", required=False)

    args = parser.parse_args()
    
    #Function Exec calls
    if args.command == 'simdata':
        vcf_to_hap(args.i,args.hpath,args.hapout, args.phenout)
    elif args.command == 'gwas':
        readData(args.geno,args.pheno,args.out,args.maf)
    else:
        parser.print_help()
        
main()
