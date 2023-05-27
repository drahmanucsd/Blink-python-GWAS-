import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

# Read the GWAS dataset into a dataframe
df = pd.read_csv('fake_gwas_data2.csv')

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
