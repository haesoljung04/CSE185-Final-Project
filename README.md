# myGWAS: A Genome-Wide Association Study Tool

## Description

myGWAS is a command-line tool designed to perform GWAS on given VCF and phenotype files, focusing on linear regression analysis and providing visualization of the results through Manhattan and QQ plots.

## Installation Instructions
NOTE: currently only working on datahub so please test on there...

First and foremost, just `git clone` this whole thing.\
To install the required libraries contained in the `requirements.txt` file, you can use `pip`. For example:

```bash
pip install -r requirements.txt
```

To install `myGWAS`, use the following command:

```bash
python setup.py install
```

If you lack root access, use the `--user` option:

```bash
pip install --user -r requirements.txt
python setup.py install --user
```

After installation, you can check the usage instructions by running(make sure you are in the myGWAS directory):

```bash
chmod +x myGWAS.py
./myGWAS.py --help
```

## Basic Usage

The basic format for running `myGWAS` is:

```bash
python myGWAS.py --vcf <path_to_vcf_file> --pheno <path_to_phenotype_file> --out <output_file_prefix> --linear --maf 0.05 --allow-no-sex
```

Example commands using sample files in the `example-files` directory(make sure you are in the CSE185-Final-Project directory when running these commands):

```bash
python ~/CSE185-Final-Project/myGWAS/myGWAS.py --vcf ~/CSE185-Final-Project/example-files/example1.vcf --pheno ~/CSE185-Final-Project/example-files/example1.phen --out example1 --linear --maf 0.05 --allow-no-sex
```

```bash
python ~/CSE185-Final-Project/myGWAS/myGWAS.py --vcf ~/CSE185-Final-Project/example-files/pruned_coatColor_maf_geno.vcf.gz --pheno ~/CSE185-Final-Project/example-files/coatColor.phen --out example1 --linear --maf 0.05 --allow-no-sex
```

## Complete Usage Instructions

### Required Arguments
The required input to the GWAS command is a VCF file and a Phenotype file. Users could also specify additional options as listed below:

- `--vcf FILE` : The inputted VCF file contains genetic variant data.
- `--pheno FILE` : The phenotye file contains trait information for the samples.
  
Option Parameters:
- `-o FILE, --out FILE` : This command writes the output to a file. By default, the output is typically written to stdout.
- `--linear` : This command performs linear regression analysis, specifying the type of statistical test that will be used in the GWAS command
- `--maf FLOAT` : This command filters the variants by the MAF (minor allele frequency). Only certain variants with float value greater than the threshold will be incorporated into the analysis, meaning that variants with MAF below the threshold will be excluded.
- `-- allow-no-sex` : This command will permit the analysis to run smoothly 

### Example Commands

```bash
python ~/CSE185-Final-Project/myGWAS/myGWAS.py --vcf ~/CSE185-Final-Project/example-files/pruned_coatColor_maf_geno.vcf.gz --pheno ~/CSE185-Final-Project/example-files/coatColor.phen --out example1 --linear --maf 0.05 --allow-no-sex
```

## File Format

The output file format for my `GWAS` tool is equivalent to the standard GWAS method. For a more detail specification, see 

## Contributors
This repository was written by Lydia Roh (A17002778), Haesol Jung (A17348180), and Meiqi Lai (A17043227).

Please submit a pull request in the event that you have any corrections or suggestions. Thank you!
