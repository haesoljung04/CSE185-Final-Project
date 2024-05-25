# myGWAS: A Genome-Wide Association Study Tool

## Description

myGWAS is a command-line tool designed to perform GWAS on given VCF and phenotype files, focusing on linear regression analysis and providing visualization of the results through Manhattan and QQ plots.

## Installation Instructions
NOTE: currently only working on datahub so please test on there...
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


### Example Commands


## File Format


## Contributors
Lydia Roh (A17002778), Haesol Jung (A17348180), Meiqi Lai (A17043227)
