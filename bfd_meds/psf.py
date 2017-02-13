#Standard library
import re
import glob
import collections

#Internal

#External
import galsim
import galsim.des
import astropy.io.fits

class PsfSource(object):
    pass

class MissingPSFError(StandardError):
    pass

class TooManyPSFsError(StandardError):
    pass



class LimitedSizeDict(collections.OrderedDict):
    """
    A dictionary with a maximum size, suitable as a simple
    Create as with a normal OrderedDict with the additional size_limit keyword arg.

    From http://stackoverflow.com/questions/2437617/limiting-the-size-of-a-python-dictionary
    """
    def __init__(self, *args, **kwds):
        self.size_limit = kwds.pop("size_limit", None)
        collections.OrderedDict.__init__(self, *args, **kwds)
        self._check_size_limit()

    def __setitem__(self, key, value):
        collections.OrderedDict.__setitem__(self, key, value)
        self._check_size_limit()

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)



class PsfexSource(PsfSource):
    def __init__(self, cache_size):
        self.cache = LimitedSizeDict(size_limit=cache_size)

    def _get_psf_data(self, fits_filename, hdu_name):
        data = self.cache.get(hdu_name)

        if data is None:
            fitsfile = astropy.io.fits.open(fits_filename)
            hdu = fitsfile[hdu_name]
            data = galsim.des.DES_PSFEx(hdu)

        self.cache[hdu_name] = data

        return data


    def evaluate_psfex(self, psfex, x_image, y_image, nsidex=32, nsidey=32, upsampling=1, offset=None):
        """Return an image of the PSFEx model of the PSF as a np array.
        Stolen from Barney Rowe many years ago

        Arguments
        ---------
        psfex       A galsim.des.PSFEx instance opened using, for example,
                    `psfex = galsim.des.DES_PSFEx(psfex_file_name)`.
        x_image     Floating point x position on image [pixels]
        y_image     Floating point y position on image [pixels]

        nsidex      Size of PSF image along x [pixels]
        nsidey      Size of PSF image along y [pixels]
        upsampling  Upsampling (see Zuntz et al 2013)
        distort_function  A function taking the GSObject returned by getPSF and returning a new object to be drawn.  Or None

        Returns a np array with shape (nsidey, nsidex) - note the reversal of y and x to match the
        np internal [y, x] style array ordering.  This is to ensure that `pyfits.writeto()` using the
        ouput array creates FITS-compliant output.
        """
        import galsim
        image = galsim.ImageD(nsidex, nsidey)
        #Note galsim uses 1-offset convention whereas coordinates in meds file are 0-offset:
        x_image_galsim = x_image+1
        y_image_galsim = y_image+1
        psf = psfex.getPSF(galsim.PositionD(x_image_galsim, y_image_galsim))

        psf.drawImage(image, scale=1.0/upsampling, offset=offset, method='no_pixel')
        return image.array


class CollectedMedsPsfexSource(PsfexSource):
    def __init__(self, template, cache_size=25):
        super(CollectedMedsPsfexSource, self).__init__(cache_size)
        self.template = template

    def _get_psfex_filename(self, tilename, band):
        return self.template.format(tilename=tilename, band=band)

    def _get_psfex_hdu(self, fitsfile, exposure, band, ccd):
        pattern = "{0}_{1}_c{2}".format(exposure, band, ccd)
        if pattern in fitsfile:
            return fitsfile[pattern]
        raise MissingPSFError(pattern)



    def get_psf(self, tilename, band, exposure, ccd, col, row, stamp_size, jacobian):
        fits_filename = self._get_psfex_filename(tilename, band)
        hdu_name = "{0}_{1}_c{2}".format(exposure, band, ccd)
        psf_data = self._get_psf_data(fits_filename, hdu_name)
        psf_image = self.evaluate_psfex(psf_data, col, row, stamp_size, stamp_size)
        return psf_image



class DirectoryPsfexSource(PsfexSource):
    def __init__(self, directory, cache_size=25):
        super(DirectoryPsfexSource, self).__init__(cache_size)
        self.directory = directory

    def _get_psfex_filename(self, exposure, band, ccd):
        pattern = "{0}/{1}_{2}_c{3}_*_psfexcat.psf".format(self.directory, exposure, band, ccd)
        files = glob.glob(pattern)
        if len(files)==0:
            raise MissingPSFError(pattern)
        elif len(files)>1:
            raise TooManyPSFsError(pattern)
        return files[0]

    def get_psf(self, tilename, band, exposure, ccd, col, row, stamp_size, jacobian):
        fits_filename = self._get_psfex_filename(exposure, band, ccd)
        hdu_name = "PSF_DATA"
        psf_data = self._get_psf_data(fits_filename, hdu_name)
        psf_image = self.evaluate_psfex(psf_data, col, row, stamp_size, stamp_size)
        return psf_image




