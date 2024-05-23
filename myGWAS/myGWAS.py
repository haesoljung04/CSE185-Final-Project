#!/opt/homebrew/bin/python3

import argparse

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

if __name__ == "__main__":
    main()
