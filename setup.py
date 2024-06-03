from setuptools import setup, find_packages


setup(
  name='myGWAS',
  version='0.1.0',
  description='CSE185 Final Project',
  author='Haesol Jung, Meiqi Lai, Lydia Roh',
  author_email='haj008@ucsd.edu, melai@ucsd.edu, lroh@ucsd.edu',
  packages=find_packages(),
  install_requires=[
        "cyvcf2",
        "scipy",
        "pandas",
        "seaborn",
        "matplotlib>=3.5.0",
        "adjustText",
        "setuptools", 
  ],
  entry_points={
    "console_scripts": [
        "myGWAS=myGWAS.myGWAS:main"
    ],
  },
)
  
