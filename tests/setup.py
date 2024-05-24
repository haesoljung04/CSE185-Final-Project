setup(
  name='myGWAS',
  version=VERSION,
  description='CSE185 Final Project',
  author='Haesol Jung, Meiqi Lai, Lydia Roh',
  author_email='haj008@ucsd.edu, melai@ucsd.edu, lroh@ucsd.edu',
  packages=find_packages(),
  entry_points={
    "console_scripts": [
        "myGWAS=myGWAS.myGWAS:main"
    ],
  },
)
  
