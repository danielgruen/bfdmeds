import ngmixer
import meds
import fitsio

# this is an example -- download tile first into data/

nstamps=100

lf=["test_meds_subtracted_replace_bad_reject_outliers.fits", "test_meds_subtracted_reject_outliers.fits", "test_meds_reject_outliers.fits"]

# (1) extract, in three different ways

m=ngmixer.imageio.extractor_corrector.MEDSExtractorCorrector("data/DES2356-5705_r2590p01-test-mof-002.fits", "data/DES2356-5705_r2590p01-r-meds-Y3A1.fits.fz", 0, nstamps-1, lf[0], reject_outliers=True, replace_bad=True)

m=ngmixer.imageio.extractor_corrector.MEDSExtractorCorrector("data/DES2356-5705_r2590p01-test-mof-002.fits", "data/DES2356-5705_r2590p01-r-meds-Y3A1.fits.fz", 0, nstamps-1, lf[1], reject_outliers=True, replace_bad=False)

m=meds.MEDSExtractor("data/DES2356-5705_r2590p01-r-meds-Y3A1.fits.fz", 0, nstamps-1, lf[2])

# (2) write out a mosaic image so we can inspect those


for f in lf:
  m=meds.MEDS(f)
  for o in range(nstamps):
    for t in ["image","weight","seg"]:
      data=m.get_mosaic(o, type=t)
      fitsio.write(t+"_"+str(o)+"_"+f,data)
