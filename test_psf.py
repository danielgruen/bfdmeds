import sys
import os
# On JZ's machine:  python test_psf.py /Users/jaz/data/DES2348-5831/DES2348-5831_r2590p01_r_meds-Y3A1.fits.fz

if len(sys.argv)!=2:
	print "Run: python test_psf.py /path/to/medsfile"
	print "Must be PSFEx files in the same directory"
	sys.exit(1)

filename = sys.argv[1]

import bfd_meds

import matplotlib
matplotlib.use("agg")
import pylab


dirname, _ = os.path.split(filename)

# Set up the sources of PSFs - in this case a directory of files
psf_source = bfd_meds.DirectoryPsfexSource("/Users/jaz/data/DES2348-5831")
# Open the MEDS file itself
m = bfd_meds.BFDMEDS(filename, psf_source)

# Get and save a picture of a PSF
psf = m.get_psf(20, 1)
pylab.imshow(psf, interpolation='nearest')
pylab.savefig("psf.png")