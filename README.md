# CSE185-Final-Project
cs185sp24 final project. Collaborators: Meiqi Lai, Lydia Roh

## Installation Instructions
Installation requires all libraries contained within the `requirements.txt` file to be installed. All of this can be accomplished through `pip`.
For example: `pip install cyvcf2`

Once you have installed all required libraries, proceed to install `myGWAS`:
`python setup.py install`

If you lack root access, please use the options specified in the following commands instead:
`pip install --user cyvcf2`
`python setup.py install --user `

After all this, use `myGWAS --help` to learn how to use this tool!
## Basic Usage
The basic format of `myGWAS` is:
`myGWAS --vcf example_gwas.vcf --pheno example_gwas.phen --out example_gwas --linear --maf 0.05 --allow-no-sex`

To run `myGWAS` on a small test example(using files in the example-files directory on this repo):
`python3 myGWAS.py --vcf example1.vcf --pheno example1.phen --out example1 --linear --maf 0.05 --allow-no-sex`

