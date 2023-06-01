from haptools import data
from haptools.sim_phenotype import Haplotype
import argparse
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import os

def vcf_to_hap(str vcf_path):
    # which variants do we want to write to the haplotype file?
    variants = {"rs149635655", "rs141306699"}
    #df = pd.read_csv("/home/jtliew/teams/31/variants.csv", header=None, names=['Variants'])
    #df = pd.read_csv("/home/jtliew/teams/31/variants.csv")
    #val = df.get("Variants")

    #for i in val:
        variants.add(i)

    #print(variants)
    
    # load the genotypes file
    # you can use either a VCF or PGEN file
    gt = data.GenotypesVCF(vcf_path)
    gt.read(variants=variants)

    # initialize an empty haplotype file
    hp = data.Haplotypes("output.hap", haplotype=Haplotype)
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
    
    
def linear_reg(str data_path):
    # Read the GWAS dataset into a dataframe
    #df = pd.read_csv('fake_gwas_data2.csv')
    df = pd.read_csv(data_path)

    # Extract the genotypes and phenotype
    genotypes = df[['Genotype1', 'Genotype2', 'Genotype3']]
    phenotype = df['Phenotype']

    # Encode the genotypes as numeric values
    genotypes_encoded = genotypes.replace({'A/A': 0, 'A/T': 1, 'A/C': 2, 'G/T': 3,
                                           'T/T': 4, 'T/C': 5, 'T/G': 6,
                                           'C/C': 7, 'C/G': 8, 'G/G': 9})

    # Add an intercept term to the genotypes
    genotypes_encoded = sm.add_constant(genotypes_encoded)

    # Fit a linear regression model
    model = sm.OLS(phenotype, genotypes_encoded)
    results = model.fit()

    # Perform GWAS
    p_values = results.pvalues[1:]  # Exclude the intercept term
    significant_snps = p_values[p_values < 0.05]  # Adjust the significance threshold as needed

    # Print significant SNPs
    for snp, p_value in significant_snps.items():
        print(f'SNP: {snp}, p-value: {p_value}')

    print(results.summary())

    # Predict the phenotype values
    predicted_phenotype = results.predict(genotypes_encoded)

    # Create a scatter plot of actual vs predicted phenotype values
    plt.scatter(phenotype, predicted_phenotype)
    plt.xlabel('Actual Phenotype')
    plt.ylabel('Predicted Phenotype')
    plt.title('Linear Regression - Actual vs Predicted Phenotype')
    plt.show()

def main():
    parser = argparse.ArgumentParser(prog='blink',description="Command-line tool to perform basic GWAS")

    parser.add_argument('genotypes',type=str,help="Path to the genotype file (.vcf.gz format)")
    
    parser.add_argument('test_csv',type=str,help="Path to the test file (.csv format)")
    
    parser.add_argument('output',type=str,help="Enter your first name")

    args = parser.parse_args()
    
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand')

    vcf_to_hap = subparsers.add_parser('vcf_to_hap', help='Runs vcf_to_hap')
    linear_reg = subparsers.add_parser('linear_reg', help='Run linear_reg')

    args = parser.parse_args()

    if args.subcommand == 'vcf_to_hap':
        vcf_to_hap(args.genotypes)
    elif args.subcommand == 'linear_reg':
        linear_reg(args.test_csv)
    else:
        parser.print_help()

if __name__ == '__main--':
    main()