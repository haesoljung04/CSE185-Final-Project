#!/usr/bin/env python

# import for parsing command line arguments
import argparse

# imports for linear regression analysis
from cyvcf2 import VCF
from scipy.stats import linregress
from tqdm import tqdm

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
      print("Finished Linear Regression")
      plot(output_file + ".assoc.linear", 'manhattan_plot.png', 'qq_plot.png')
      print("Finished Plotting")
      if print_summary:
        summary(output_file + ".assoc.linear")
        print("Finished Summary")
  else:
      print("Currently, only linear regression is implemented. Use --linear please.")

# This function creates an index file to speed up linear regression
def create_index_file(vcf_file, index_file):
    # First parse VCF file using cyvcf2
    vcf = VCF(vcf_file)
    samples = vcf.samples
    
    # Create a DataFrame to store sample IDs and indices
    index_df = pd.DataFrame({"IndID": samples, "Index": range(len(samples))})
    
    # Save DataFrame to index file
    index_df.to_csv(index_file, sep="\t", header=False, index=False)

# Perform linear regression analysis(the heart of GWAS
def linear_regression(vcf_file, pheno_file, output_file, maf_threshold, allow_no_sex):
    # Parse VCF file using cyvcf2
    vcf = VCF(vcf_file)
    samples = vcf.samples
    
    # Count the total number of variants
    total_variant_count = sum(1 for _ in VCF(vcf_file))
    
    # Make Index File
    create_index_file(vcf_file, "index_file.txt")
    # Read index file
    index_df = pd.read_csv("index_file.txt", delim_whitespace=True, header=None, names=["IndID", "Index"])
    # Create a dictionary mapping sample IDs to their corresponding indices
    index_dict = dict(zip(index_df["IndID"], index_df["Index"]))

    # Read phenotype file
    pheno_df = pd.read_csv(pheno_file, delim_whitespace=True, header=None, names=["famID", "IndID", "Phenotype"])
    # Create a dictionary where each IID is a key and the value is pheno value
    pheno_dict = dict(zip(pheno_df["IndID"], pheno_df["Phenotype"]))
    unique_phenotypes = pheno_df["Phenotype"].unique()
    # Create a binary mapping based on the unique phenotypes
    binary_mapping = {phenotype: index for index, phenotype in enumerate(unique_phenotypes)}

    # Prepare output file
    output = open(output_file + ".assoc.linear", "w")
    output.write("CHR\tSNP\tBP\tA1\tTEST\tNMISS\tBETA\tSTAT\tP\n") # this will be the header

    # Initialize tqdm progress bar with total variant count
    with tqdm(total=total_variant_count, desc="Performing Linear Regression") as pbar:
        # Iterate through each variant using cyvcf2
        for variant in vcf:
            # Update progress bar
            pbar.update(1)
            # Calculate the maf
            allele_counts = variant.gt_alt_freqs
            for allele_count in allele_counts:
                maf = np.min(allele_count) / np.sum(allele_count)
                # If the maf is too small then do not include
                if maf < maf_threshold:
                    continue
                
            # Prep data for linear regression
            genotype_data = []
            phenotype_data = []
            for sample in samples:
                # Check if the sample actually exists in the pt dictionary
                if sample in pheno_dict:  
                    index = index_dict[sample]
                    gt = variant.genotypes[index]
                    genotype_data.append(sum(gt))
                    phenotype_data.append(binary_mapping[pheno_dict[sample]])
            # If we have samples and everything was correctly initialized
            if len(genotype_data) > 0:
              # Convert to np arrays for speed
              genotype_data = np.array(genotype_data, dtype=int)
              phenotype_data = np.array(phenotype_data, dtype=int)
              # Check the shape of the arrays
        
              # Do linear regression using scipy
              slope, intercept, r_value, p_value, std_err = linregress(genotype_data, phenotype_data)
              # Write results to output(formatting to 4 decimal places)
              output.write(f"{variant.CHROM}\t{variant.ID}\t{variant.POS}\t{variant.REF}\tADD\t{len(genotype_data)}\t{slope:.4f}\t{r_value / std_err:.4f}\t{p_value:.4g}\n")

    output.close()


# Plotting the Manhattan and QQ plots
def plot(linearFile, manhattan_file, qq_file):
    # Load your GWAS results
    df = pd.read_csv(linearFile, delim_whitespace=True)
  
    # Set font properties (if needed)
    # plt.rcParams['font.family'] = 'Arial'
    # plt.rcParams['font.size'] = 12
    # plt.rcParams['font.weight'] = 'normal'
  
    # sns.set(font='Arial', style='white')
  
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
        ax.set_xticklabels(chrom_df.index, fontsize=12)
  
        # Annotate the top 5 most significant SNPs
        for idx, row in top_5_snps.iterrows():
            texts.append(ax.text(row['i'], row['-log10(P)'], row['SNP'], fontsize=10, ha='right', va='bottom'))
  
    # Adjust the text positions to minimize overlap
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray'))
  
    # Set title
    plot.fig.suptitle('Manhattan plot', fontsize=16, weight='bold')
  
    # Adjust layout to remove white space
    plt.tight_layout()
  
    # Save the Manhattan plot to a file
    with tqdm(desc="Saving Manhattan plot") as pbar:
        plt.show()
        plot.savefig(manhattan_file)
        pbar.update(1)
    
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
  
    # Save the QQ plot to a file
    with tqdm(desc="Saving QQ plot") as pbar:
        plt.show()
        fig.savefig(qq_file)
        pbar.update(1)


# Output the # genome wide significant SNPS
def summary(linearFile):
    sig_SNP_count = 0
    # Open the linear file for reading
    with open(linearFile, 'r') as file:
        # Skip the header line
        next(file)
        
        # Read each line in the file
        for line in tqdm(file, desc="Calculating summary statistics"):
            fields = line.strip().split('\t')
            p_value = float(fields[8])
            # Check if the SNP is genome-wide significant
            if p_value < 5e-8:
                sig_SNP_count += 1
    
    # Print total number of significant SNPs found
    print(f"Total number of significant SNPs found: {sig_SNP_count}")



if __name__ == "__main__":
    main()
