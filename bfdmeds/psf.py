#Standard library
import re
import glob
import collections

#Internal

#External
import galsim
import galsim.des
import astropy.io.fits
try:
    import piff
except:
    print("piff not found, must install if using piff psfs")
class PsfSource(object):
    """
    A Base class representing sources of PSF information.

    Subclasses are expected to override PsfSource.get_psf
    See their docstrings for more info.
    """
    def get_psf(self, tilename, band, exposure, ccd, col, row, stamp_size, jacobian):
        raise NotImplementedError("Need to implement get_psf in PsfSource objects")


class MissingPSFError(StandardError):
    """Error used when no PSF could be found for an exposure"""
    pass

class TooManyPSFsError(StandardError):
    """Error used when more than one PSF could be found for an exposure"""
    pass



class LimitedSizeDict(collections.OrderedDict):
    """
    A dictionary with a maximum size, suitable as a simple cache
    
    Parameters
    ----------
    Usual dictionary creation parameters plus:

    psf_size: Maximum size after which objects are removed

    From http://stackoverflow.com/questions/2437617/limiting-the-size-of-a-python-dictionary
    """
    def __init__(self, *args, **kwds):
        self.size_limit = kwds.pop("size_limit", None)
        collections.OrderedDict.__init__(self, *args, **kwds)
        self._check_size_limit()

    def __setitem__(self, key, value):
        #Remove and re-insert to update the cache
        self.pop(key,None)
        collections.OrderedDict.__setitem__(self, key, value)
        self._check_size_limit()

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)



class PsfInfoSource(PsfSource):
    """
    Base class for Piff or PSFEx-related sources of PSF information.

    Defines one method for internal use, _get_psf_data, which opens
    a PSFEx FITS file and gets and caches PSF information from it.

    Parameters
    ----------
    cache_size: number of instances of galsim.des.DES_PSFEx cached.

    """
    def __init__(self, cache_size):
        self.cache = LimitedSizeDict(size_limit=cache_size)

    def _get_psf_data(self, fits_filename, hdu_name):
        """
        Check for the named hdu in the cache.  If found, return the 
        cached HDU data.

        Otherwise open the named fits file, get the data from the named
        extension, and turn into a galsim.des.DES_PSFEx object.
        """
        data = self.cache.get((fits_filename,hdu_name))

        if data is None:
            if self.psf_type == 'psfex':
                fitsfile = astropy.io.fits.open(fits_filename)
                hdu = fitsfile[hdu_name]
                data = galsim.des.DES_PSFEx(hdu)

            if self.psf_type == 'piff':
                data = piff.read(fits_filename)

        self.cache[(fits_filename,hdu_name)] = data

        return data


    def evaluate_psf(self, psfdata, x_image, y_image, nsidex=32, nsidey=32, upsampling=1, jacobian=None, offset=None, return_psf_object=False):
        """Return an image of the Piff or PSFEx model of the PSF as a np array.
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
        wcs=None
        if jacobian is not None:
            wcs=galsim.AffineTransform(jacobian['dudcol'],jacobian['dudrow'],jacobian['dvdcol'],jacobian['dvdrow'],origin=galsim.PositionD(jacobian['col0'],jacobian['row0']),world_origin=galsim.PositionD(0,0))

        # setup image with desired wcs of postage stamp
        image = galsim.ImageD(nsidex, nsidey,wcs=wcs)

        #Note galsim uses 1-offset convention whereas coordinates in meds file are 0-offset:
        x_image_galsim = x_image+1
        y_image_galsim = y_image+1

        if self.psf_type == 'psfex':
            # draw psf in correct WCS to image
            psf = psfdata.getPSF(galsim.PositionD(x_image_galsim, y_image_galsim))
            psf_world=wcs.toWorld(psf)
            psf_world.drawImage(image, offset=offset, method='no_pixel')

        if self.psf_type == 'piff':
            # draw psf to image with desired WCS
            image=psfdata.draw(x=x_image_galsim, y=y_image_galsim, image=image)
            psf_world=psfdata

        if return_psf_object:
            return (image, psf_world)
        else:
            return image


'''
class CollectedMedsPsfexSource(PsfexSource):
    """
    Class to represent PSFEx files collated into a single file
    as a number of extensions.

    Implemented get_psf so can be used by BFDMEDS.

    Parameters
    ----------
    template: a string to be used as a template for filenames
        for a given tilename and band, e.g. /path/to/psfs/{tilename}_{band}.fits

    """
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



    def get_psf(self, tilename, band, exposure, ccd, col, row, stamp_size, jacobian,return_image=True):
        """
        Get an image of the PSF for the given exposure and position

        Parameters
        ----------
        tilename:
            string, Identifier for the tile, used to find the correct collected PSF file
        band:
            string, Identifier for the band, used to find correct PSF
        exposure:
            string, Identifier for the exposure
        ccd:
            string or int, Identifier for the chip
        col:
            integer in MEDS indexing of the chip column to evaluate at
        row:
            integer in MEDS indexing of the chip row to evaluate at
        stamp_size:
            integer size of the psf image to return in pixels. Since PSFEx models
            are stored in pixel space this implies a relative resolution
        jacobian:
            unused for now

        Returns:
            A numpy stamp_size x stamp_size image of the PSF

        """
        fits_filename = self._get_psfex_filename(tilename, band)
        hdu_name = "{0}_{1}_c{2}".format(exposure, band, ccd)
        psf_data = self._get_psf_data(fits_filename, hdu_name)
        if return_image:
 #           psf_image = self.evaluate_psfex(psf_data, col, row, stamp_size, stamp_size, offset=(0.5,0.5),return_image=True)
            psf_image = self.evaluate_psfex(psf_data, col, row, stamp_size, stamp_size,return_image=True)
            return psf_image
        else:
            psf_obj = self.evaluate_psfex(psf_data, col, row, stamp_size, stamp_size,return_image=False)
            return psf_obj
        # offset so central pixel is (stamp_size/2,stamp_size/2) when starting at 0
'''



class DirectoryPsfInfoSource(PsfInfoSource):
    """
    Class to represent Piff or PSFEx files stored in a single directory.
    Implemented get_psf so can be used by BFDMEDS.

    This version currently assumes it is isntantiated afresh for 
    each tilename, so that the tile name is not used to form 
    the path to search.

    Parameters
    ----------
    directory: The name of the directory to search for .psf files.

    """
    def __init__(self, directory, psf_type, cache_size=25):

        super(DirectoryPsfInfoSource, self).__init__(cache_size)
        self.directory = directory
        self.psf_type = psf_type

    def _get_psf_filename(self, exposure, band, ccd):

        if self.psf_type == 'psfex':
            pattern = "{0}/{1}_{2}_c{3}_*_psfexcat.psf".format(self.directory, exposure, band, ccd)
        elif self.psf_type == 'piff':
            pattern = "{0}/{1}_{2}_c{3}_*_piff.fits".format(self.directory, exposure, band, ccd)
        else:
            raise Exception("incorrect psf_type given (not piff or psfex)")

        files = glob.glob(pattern)
        if len(files)==0:
            raise MissingPSFError(pattern)
        elif len(files)>1:
            raise TooManyPSFsError(pattern)
        return files[0]


    def _get_psfex_coadd_filename(self, tilename, ccd,  request_attempt, band):
        pattern = "{0}/{1}_{2}_{3}_{4}_psfcat.psf".format(self.directory, tilename, ccd, request_attempt, band)
        files = glob.glob(pattern)
        if len(files)==0:
            raise MissingPSFError(pattern)
        elif len(files)>1:
            raise TooManyPSFsError(pattern)
        return files[0]

    def get_psf(self, tilename, band, exposure, ccd, col, row, stamp_size, jacobian,return_psf_object=False):
        """
        Get an image of the PSF for the given exposure and position

        Parameters
        ----------
        tilename:
            string, Identifier for the tile, currently unused
        band:
            string, Identifier for the band, used to find correct PSF
        exposure:
            string, Identifier for the exposure
        ccd:
            string or int, Identifier for the chip
        col:
            integer in MEDS indexing of the chip column to evaluate at
        row:
            integer in MEDS indexing of the chip row to evaluate at
        stamp_size:
            integer size of the psf image to return in pixels. Since PSFEx models
            are stored in pixel space this implies a relative resolution
        jacobian:
            unused for now

        Returns:
            A numpy stamp_size x stamp_size image of the PSF

        """
        fits_filename = self._get_psf_filename(exposure, band, ccd)
        hdu_name = "PSF_DATA"
        psf_data = self._get_psf_data(fits_filename, hdu_name)

        psf_image = self.evaluate_psf(psf_data, col, row, nsidex=stamp_size, nsidey=stamp_size,jacobian=jacobian,return_psf_object=return_psf_object) 

        return psf_image



    def get_coadd_psf(self, tilename, band, ccd, request_attempt, col, row, stamp_size, jacobian,return_psf_object=False):
        """
        Get an image of the PSF for the given exposure and position

        Parameters
        ----------
        tilename:
            string, Identifier for the tile, currently unused
        band:
            string, Identifier for the band, used to find correct PSF
        exposure:
            string, Identifier for the exposure
        ccd:
            string or int, Identifier for the chip
        col:
            integer in MEDS indexing of the chip column to evaluate at
        row:
            integer in MEDS indexing of the chip row to evaluate at
        stamp_size:
            integer size of the psf image to return in pixels. Since PSFEx models
            are stored in pixel space this implies a relative resolution
        jacobian:
            unused for now

        Returns:
            A numpy stamp_size x stamp_size image of the PSF

        """
        fits_filename = self._get_psfex_coadd_filename(tilename, ccd,  request_attempt, band)
        hdu_name = "PSF_DATA"
        psf_data = self._get_psf_data(fits_filename, hdu_name)

        psf_image = self.evaluate_psf(psf_data, col, row, nsidex=stamp_size, nsidey=stamp_size, jacobian=jacobian,return_psf_object=return_psf_object) 
        return psf_image




