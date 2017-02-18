import os
import glob
from distutils.core import setup, Extension
import numpy

dependencies = ['numpy', 'meds', 'astropy', 'pixmappy']

setup(name="bfdmeds", 
      version="0.1",
      description="Wrapper around MEDS files for interfacing to BFD code and new astrometry",
      license = "GNU GPLv3",
      author="Daniel Gruen",
      author_email="dgruen@stanford.edu",
      url="https://github.com/danielgruen/bfdmeds",
      packages=['bfdmeds'],
      install_requires=dependencies)


