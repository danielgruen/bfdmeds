import bfdmeds

# meds=bfdmeds.BFDMEDS("data/DES2348-5831_r2590p01_r_meds-Y3A1.fits.fz")
## works, but you won't be able to get stamps with neighbors subtracted

meds=bfdmeds.BFDMEDS("data_example/test_meds.fits", "data_example/test_meds_subtracted.fits", "data_example/test_meds_subtracted_replace_bad.fits")
# give three MEDS files, regular, neighbor subtracted, neighbor subtracted + bad-pixel filled

print meds.get_cutout_list(0, skip_coadd=True)
# all get_*_list methods take an optional skip_coadd parameter; if true, the coadd stamp / coordinates / etc. are not in the returned list
