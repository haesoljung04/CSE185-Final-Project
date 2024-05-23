#!/opt/homebrew/bin/python3

import argparse

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

def plot():
  # Load your GWAS results
  df = pd.read_csv("lab3_gwas.assoc.linear", delim_whitespace=True)
  
  # Set font properties
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



if __name__ == "__main__":
    main()
