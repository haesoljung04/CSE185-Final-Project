#!/usr/bin/env python

# import for parsing command line arguments
import argparse

# imports for linear regression analysis
from cyvcf2 import VCF
from scipy.stats import linregress

# imports for manhattan plot/qq plot
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text

def main():
  parser = argparse.ArgumentParser(
      prog="myGWAS",
      description="Commmand line script to perform GWAS of given VCF and Phenotype Files"
  )

  # Required Arguments
  parser.add_argument('--vcf', type=str, required=True, help='Path to VCF file')
  parser.add_argument('--pheno', type=str, required=True, help='Path to phenotype file')
  
  # Optional arguments
  parser.add_argument('-o', '--out', type=str, help='Write output to File (Default: stdout)',\
              required=False, metavar='FILE')
  parser.add_argument('--linear', action='store_true', help='Perform linear regression',\
              required=False)
  parser.add_argument('--maf', type=float, default=0.05, help='Minor allele frequency threshold (default: 0.05)',\
              required=False)
  parser.add_argument('--allow-no-sex', action='store_true', help='Allow samples without sex information',\
              required=False)
  parser.add_argument('--summary', action='store_true', help='Print summary of statistics to a file')
  
  args = parser.parse_args()

  # Accessing all the command line arguments
  vcf_file = args.vcf
  pheno_file = args.pheno
  output_file = args.out if args.out else "myGWAS"  # the default output file will be named myGWAS.assoc.linear
  perform_linear = args.linear
  maf_threshold = args.maf
  allow_no_sex = args.allow_no_sex
  print_summary = args.summary
  if perform_linear:
      linear_regression(vcf_file, pheno_file, output_file, maf_threshold, allow_no_sex)
      plot(output_file + ".assoc.linear")
      if print_summary:
        summary(output_file + ".assoc.linear")
  else:
      print("Currently, only linear regression is implemented. Use --linear please.")


# Perform the linear regression for quantitative traits
def linear_regression(vcf_file, pheno_file, output_file, maf_threshold, allow_no_sex):
    # Read phenotype file
    pheno_df = pd.read_csv(pheno_file, delim_whitespace=True, header=None, names=["famID", "IndID", "Phenotype"])
    # Create a dictionary where each IID is a key and the value is pheno value
    pheno_dict = dict(zip(pheno_df["IndID"], pheno_df["Phenotype"]))
    # Parse VCF file using cyvcf2
    vcf = VCF(vcf_file)
    samples = vcf.samples

    # Prepare output file
    output = open(output_file + ".assoc.linear", "w")
    output.write("CHR\tSNP\tBP\tA1\tTEST\tNMISS\tBETA\tSTAT\tP\n") # this will be the header
    # Iterate through each variant using cyvcf2
    for variant in vcf:
        # Calculate MAF to check if we include the variant in our calculations
        alleles = variant.gt_alleles
        allele_counts = np.zeros(2, dtype=int)
        allele_counts = np.sum(gt_bases.astype(int), axis=0)
        maf = np.min(allele_counts) / np.sum(allele_counts)
        # If the maf is too small then do not include
        if maf < maf_threshold:
            continue
        
        # Prep data for linear regression
        genotype_data = []
        phenotype_data = []
        for i, sample in enumerate(samples):
            # Check if the sample actually exists in the phenotype dictionary
            if sample in pheno_dict:
                gt = variant.genotypes[i][0:2]
                genotype_data.append(sum(gt))
                phenotype_data.append(pheno_dict[sample])
        # If we have samples and everything was correctly initialized
        if len(genotype_data) > 0:
            genotype_data = np.array(genotype_data)
            phenotype_data = np.array(phenotype_data)

            # Do linear regression using scipy
            slope, intercept, r_value, p_value, std_err = linregress(genotype_data, phenotype_data)
            
            # Write results to output(formatting to 4 decimal places)
            output.write(f"{variant.CHROM}\t{variant.ID}\t{variant.POS}\t{variant.REF}\tADD\t{len(genotype_data)}\t{slope:.4f}\t{r_value / std_err:.4f}\t{p_value:.4g}\n")

    output.close()

  

# Plotting the Manhattan and QQ plots
def plot(linearFile):
  # Load your GWAS results
  df = pd.read_csv("linearFile", delim_whitespace=True)
  
  # Set font properties(not working rn)
  plt.rcParams['font.family'] = 'Arial'
  plt.rcParams['font.size'] = 12
  plt.rcParams['font.weight'] = 'normal'
  
  sns.set(font='Arial', style='white')
  
  # Data preparation
  df['-log10(P)'] = -np.log10(df['P'])
  df['CHR'] = df['CHR'].astype('category')
  df = df.sort_values(['CHR', 'BP'])
  df.reset_index(inplace=True, drop=True)
  df['i'] = df.index

  # Coloring
  unique_chromosomes = sorted(df['CHR'].unique())
  num_chromosomes = len(unique_chromosomes)
  colors = sns.color_palette("husl", num_chromosomes)
  
  # Generate Manhattan plot
  plot = sns.relplot(data=df, x='i', y='-log10(P)', aspect=3.7, 
                     hue='CHR', palette=colors, legend=None, marker='o', alpha=1, s=15)  # Adjust marker type and transparency
  
  # Identify the top 5 most significant SNPs
  top_5_snps = df.nsmallest(5, 'P')
  
  # Iterate over each Axes object in the FacetGrid
  texts = []
  for ax in plot.axes.flat:
      # Add a horizontal line for genome-wide significance
      ax.axhline(y=-np.log10(5e-8), color='black', linestyle='--', linewidth=1)
      
      # Set x-axis label and tick labels
      chrom_df = df.groupby('CHR')['i'].median()
      ax.set_xlabel('Chromosome')
      ax.set_xticks(chrom_df)
      ax.set_xticklabels(chrom_df.index, fontname='Georgia', fontsize=12)
  
      # Annotate the top 5 most significant SNPs
      for idx, row in top_5_snps.iterrows():
          texts.append(ax.text(row['i'], row['-log10(P)'], row['SNP'], fontsize=10, ha='right', va='bottom', fontname='Georgia'))
  
  # Adjust the text positions to minimize overlap
  adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray'))
  
  # Set title
  plot.fig.suptitle('Manhattan plot', fontsize=16, weight='bold', fontname='Georgia')
  
  # Adjust layout to remove white space
  plt.tight_layout()
  
  # Show the Manhattan plot
  plt.show()
  
  # Generate the QQ plot
  fig, ax = plt.subplots(figsize=(8, 8))
  
  # Sort the p-values
  sorted_p_values = np.sort(df['P'])
  # Calculate expected quantiles
  expected_quantiles = -np.log10(np.linspace(1 / len(sorted_p_values), 1, len(sorted_p_values)))
  # Plot the observed versus expected quantiles
  ax.scatter(expected_quantiles, -np.log10(sorted_p_values), color=colors[0], alpha=0.7)
  ax.plot([0, max(expected_quantiles)], [0, max(expected_quantiles)], color=colors[1], linestyle='--')
  ax.set_xlabel('Expected -log10(p-value)')
  ax.set_ylabel('Observed -log10(p-value)')
  ax.set_title('QQ Plot')
  
  plt.show()

# Output the # genome wide significant SNPS
def summary(linearFile):
    sig_SNP_count = 0
    # Open the linear file for reading
    with open(linearFile, 'r') as file:
        # Read each line in the file
        for line in file:
            fields = line.strip().split('\t')
            p_value = float(fields[8])
            # Check if the SNP is genome-wide significant
            if p_value < 5e-8:
                significant_snps_count += 1
    
    # Print total number of significant SNPs found
    print(f"Total number of significant SNPs found: {significant_snps_count}")


if __name__ == "__main__":
    main()
