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

def vcf_to_hap(vcf_path:str, hpath: str, hap_out_path: str, pheno_out_path):
#THE CODE BELOW IS CREDITED TO Massarat, A. R., Lamkin, M., Reeve, C., Williams, A. L., D’Antonio, M., & Gymrek, M. Haptools: a toolkit for admixture and haplotype analysis [Computer software]. https://github.com/CAST-genomics/haptools who created haptools (Was granted permission to use from Professor Gymrek)
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

def gwas_cmd(geno_path, pheno_path,out_path):
    print("hi")
    
def process_SNPS(genotype_df, phenotype_df, row_number):
    
    #CREATE THE ARRAY OF GENOTYPE DATA
    column_names = genotype_df.columns.tolist() #this is a list of all the column names/individuals
    values = genotype_df.iloc[row_number].values.tolist()
    #print(values)

    sum_values = [int(string.split('|')[0]) + int(string.split('|')[1]) for string in values] #this is a list of the summed values per row

    gts = np.array(sum_values)
    
    if np.var(gts) == 0:
        print(f"Zero variance for row {row_number}")
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

def plot_qq(pval_list, beta_list, genotype_data):
    zero_list = [0] * len(pval_list) #NEED TO CHANGE THIS
    #first make the dataframe for a qq plot
    qq_df = pd.DataFrame({'CHR': genotype_data.get("#CHROM"), 'SNP': genotype_data.get("ID"), 
                      "BP" : genotype_data.get("POS"), "A1" : genotype_data.get("REF"), "TEST" : genotype_data.get("FILTER"), 
                      "NMISS" : genotype_data.get("QUAL"), "BETA" : beta_list, "STAT" : zero_list, "P" : pval_list})
    
    fig, (ax0, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
    fig.set_size_inches((15, 5))
    qqman.manhattan(qq_df, ax=ax0)
    qqman.qqplot(qq_df, ax=ax1)
    plt.savefig('test.png', bbox_inches='tight') #let user choose the output file name??

def readData(genotypeData, phenotypeData):
    pval_list = [] #hold pvals
    beta_list = [] #holds beta vals
    
    lines_to_skip = 257 #change this
    #read genotype data
    genotypes = pd.read_csv(genotypeData, skiprows=lines_to_skip, sep='\t')
    #read phenotype data
    phenotypes = pd.read_csv(phenotypeData, sep='\t', header=None, names=['ID', 'ID2', 'Val']).drop(columns=["ID2"])
    #takes only the top 10 SNPS; only for testing
    #get_snps = genotypes.head(len(genotypes))
    
    #This runs our code on a small chunk of data
    get_snps = genotypes.head(100)
    get_snps_reformat = get_snps.drop(columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
    
    for i in range(get_snps_reformat.shape[0]):
        out_beta, out_pval = process_SNPS(get_snps_reformat, phenotypes, i) #CALL process_SNPS
        pval_list.append(out_pval)
        beta_list.append(out_beta)
    
    plot_qq(pval_list, beta_list, get_snps)
    #call plot_data using pval_list and beta_list

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
    #gwas.add_argument('--o',type=str,help="Path to output graph file (.png format)", required=True)

    args = parser.parse_args()
    
    #Function Exec calls
    if args.command == 'simdata':
        vcf_to_hap(args.i,args.hpath,args.hapout, args.phenout)
    elif args.command == 'gwas':
        readData(args.g,args.p)
    else:
        parser.print_help()
        
main()
