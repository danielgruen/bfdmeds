import ngmixer
import meds
import fitsio

# this is an example -- download tile first into data/

nstamps=100

# (1) extract, in three different ways

#m=ngmixer.imageio.extractor_corrector.MEDSExtractorCorrector("data/DES2356-5705_r2590p01-test-mof-002.fits", "data/DES2356-5705_r2590p01-r-meds-Y3A1.fits.fz", 1, nstamps, "test_meds_subtracted_replace_bad.fits", replace_bad=True)
# ngmixer's MEDSExtractorCorrector expects the band in the MEDS filename as, e.g., -r-; I've added _r_ in my local branch

#m=ngmixer.imageio.extractor_corrector.MEDSExtractorCorrector("data/DES2356-5705_r2590p01-test-mof-002.fits", "data/DES2356-5705_r2590p01-r-meds-Y3A1.fits.fz", 1, nstamps, "test_meds_subtracted.fits", replace_bad=False)

#m=meds.MEDSExtractor("data/DES2356-5705_r2590p01-r-meds-Y3A1.fits.fz", 1, 100, "test_meds.fits")

# (2) write out a mosaic image so we can inspect those

lf=["test_meds_subtracted_replace_bad.fits", "test_meds_subtracted.fits", "test_meds.fits"]

for f in lf:
  m=meds.MEDS(f)
  for o in range(nstamps):
    for t in ["image","weight","seg"]:
      data=m.get_mosaic(o, type=t)
      fitsio.write(t+"_"+str(o)+"_"+f,data)
