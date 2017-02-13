"""
Defines a derived class of MEDS (Multi Epoch Data Structures) to provide the
input BFD needs.

TODO: make sure Katie and I agree on an interface
TODO: implement the getters properly
TODO: figure out at which point we determine which exposures to cull (because of masking / non-white noise) and how that culling happens

"""

import meds
import numpy as np
import sys
import collections

class BFDMEDS(meds.MEDS):


    def __init__(self, filename, psf_source, filename_mofsub=None, filename_mofsub_badfill=None):
        """
        Build a BFDMEDS object from three MEDS files: the standard one, one with 
        MOF-neighbors subtracted, and one with MOF neighbors subtracted and bad pixels
        filled with central galaxy model.
        
        See extractor_corrector.py in the ngmixer repository for how to write the latter
        two versions.
        
        parameters
        ----------
        filename:
            MEDS filename
        psf_source:
            An instance of bfdpsf.PsfSource
        filename_mofsub:
            MEDS filename, MOF neighbor models subtracted from all single epoch stamps
        fiename_mofsub_badfill:
            MEDS filename, MOF neighbor models subtracted and bad pixels filled
        """
    
    
        super(BFDMEDS, self).__init__(filename)        
        
        self.psf_source = psf_source

        if(filename_mofsub is not None):
            self._meds_mofsub = meds.MEDS(filename_mofsub)
            if(self._meds_mofsub.size!=self.size):
                sys.exit("MOF subtracted and base meds file have different length. exiting.")
        else:
            self._meds_mofsub=None
        
        if(filename_mofsub_badfill is not None):
            self._meds_mofsub_badfill = meds.MEDS(filename_mofsub_badfill)
            if(self._meds_mofsub_badfill.size!=self.size):
                sys.exit("MOF subtracted bad-pixel filled and base meds file have different length. exiting.")
        else:
            self._meds_mofsub_badfill=None
            

    def get_cutout_list(self, iobj, type='image', skip_coadd=False):
         
        """
        Get an image list with all cutouts associated with this coadd object,
        subtracting neighbors and filling missing pixels in each cutout with 
        the MOF models 

        Note each individual cutout is actually a view into a larger
        mosaic of all images. The first cutout is the coadd cutout, which you 
        might want to skip for some purposes.

        parameters
        ----------
        iobj:
            Index of the object
        type: string, optional
            Cutout type. Default is 'image'. Allowed
            values are 'image','weight','seg' implemented by the base class.
            This class adds 'image_mofsub' (image with MOF neighbors subtracted)
            and 'image_mofsub_badfill' that, with the addition of bad pixels 
            being filled by model of central galaxy. 
            Types 'diff_image_mofsub' and 'diff_image_mofsub_badfill' give you the difference
            of these and the main image.
        skip_coadd:
            if True, remove the coadd postage stamp from the list (default: False)

        returns
        -------
        A list of images holding all cutouts.
        """
        
        if(skip_coadd==True):
            return self.get_cutout_list(iobj, type, False)[1:]

        if(type=='image_mofsub'):
            if(self._meds_mofsub is None):
                sys.exit("requesting MOF subtracted stamp, which was not provided")
            return self._meds_mofsub.get_cutout_list(iobj, 'image')
        if(type=='image_mofsub_badfill'):
            if(self._meds_mofsub_badfill is None):
                sys.exit("requesting MOF subtracted bad-pixel filled stamp, which was not provided")
            return self._meds_mofsub_badfill.get_cutout_list(iobj, 'image')
        if(type=='diff_image_mofsub'):
            if(self._meds_mofsub is None):
                sys.exit("requesting MOF subtracted stamp, which was not provided")
            return [a-b for a,b, in zip(super(BFDMEDS, self).get_cutout_list(iobj, 'image'),
                                        self._meds_mofsub.get_cutout_list(iobj, 'image'))]
        if(type=='diff_image_mofsub_badfill'):
            if(self._meds_mofsub_badfill is None):
                sys.exit("requesting MOF subtracted bad-pixel filled stamp, which was not provided")
            return [a-b for a,b, in zip(super(BFDMEDS, self).get_cutout_list(iobj, 'image'),
                                        self._meds_mofsub_badfill.get_cutout_list(iobj, 'image'))]
                                        
        return super(BFDMEDS, self).get_cutout_list(iobj, type)


    def get_rowcol_list(self, iobj, skip_coadd=False):
        """
        Get a list of centroid rows and columns for all cutouts associated with this
        coadd object.

        TODO: decide whether we need this since the same information is contained in
              the output of get_jacobian_list

        parameters
        ----------
        iobj:
            Index of the object
        skip_coadd:
            if True, remove the coadd rowcol from the list (default: False)

        returns
        -------
        A list of row positions of object centroid, and a list of its col positions.
        """
                
        if(skip_coadd==True):
            return self.get_rowcol_list(iobj, False)[1:]

        rowlist=[]
        collist=[]

        ncutout=self._cat['ncutout'][iobj]

        for i in xrange(ncutout):
            row,col = self.get_cutout_rowcol(iobj,i)
            rowlist.append(row)
            collist.append(col)

        return rowlist, collist


    def get_psf_list(self, psf_source, iobj, skip_coadd=False):
        """
        Get a list of PSF postage stamp images for all cutouts associated with this
        coadd object.

        parameters
        ----------
        iobj:
            Index of the object
        psf_source:

        skip_coadd:
            if True, remove the coadd PSF from the list (default: False)


        returns
        -------
        A list of PSF postage stamps.
        """

        if(skip_coadd==True):
            return self.get_psf_list(iobj, False)[1:]
            
        ncutout=self._cat['ncutout'][iobj]

        psflist=[]
        for i in xrange(ncutout):
            psflist.append(self.get_psf(iobj, i)) # TODO: implement


    def get_psf(self, iobj, icutout):
        """
        Get a a PSF image for a single exposure.

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of the exposure.  Zero for coadd.


        returns
        -------
        A numpy array containing the PSF image.
        """
        info = self.get_exposure_info(iobj, icutout)
        row, col = self.get_cutout_rowcol(iobj, icutout)
        stamp_size = self.get_cat()['box_size'][iobj]
        jacobian = None  # for now

        return self.psf_source.get_psf(info['tilename'], info['band'], info['exposure'], info['ccd'], row, col, stamp_size, jacobian)

    def get_exposure_info(self, iobj, icutout):
        """
        Get a selection of information about this exposure.

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of the exposure.  Zero for coadd.


        returns
        -------
        A numpy array containing the PSF image.
        """
        if icutout==0:
            raise ValueError("MEDS.get_exposure_info requires icutout>0 (no information available for coadds)")

        #Get the source info and thence image path
        info = self.get_source_info(iobj, icutout)
        image_path = info['image_path'] 
        # image_paths have this format:
        #"nwgint/DES2348-5831_r2590p01_D00350178_r_c03_nwgint.fits"
        # we parse it into the pieces it encodes:
        tile_part, request_attempt, exposure, band, ccd_part, _ = image_path.split("_")
        
        #Strip out the boilerplate
        ccd = ccd_part.lstrip("c")
        prefix = "nwgint/"
        tilename = tile_part[len(prefix):]

        #Return info as dictionary
        info = dict(ccd=ccd, tilename=tilename, request_attempt=request_attempt, band=band, exposure=exposure)
        return info






     
    def get_jacobian_list(self, iobj, skip_coadd=False):
        """
        Get the list of jacobians for all cutouts
        for this object.

        parameters
        ----------
        iobj:
            Index of the object
        skip_coadd:
            if True, remove the coadd jacobian from the list (default: False)
        """
        
        if(skip_coadd==True):
            return super(BFDMEDS, self).get_jacobian_list(iobj)[1:]
        
        return super(BFDMEDS, self).get_jacobian_list(iobj)


    def get_noise_list(self, iobj, skip_coadd=False):
        """
        Get a list of mean rms noise levels for all cutouts associated with this
        coadd object, which is the inverse square root of the harmonic mean of weights.
        
        parameters
        ----------
        iobj:
            Index of the object
        skip_coadd:
            if True, remove the coadd noise level from the list (default: False)
        """
        
        if(skip_coadd==True):
            return self.get_noise_list(iobj, False)[1:]

        ncutout=self._cat['ncutout'][iobj]

        wimgs=self.get_cutout_list(iobj,type='weight')

        noise=[]

        for i in xrange(ncutout):
            noise.append(1./np.sqrt(np.median(wimgs[i][wimgs[i]>0],axis=None)))

        # TODO: deal with bad/masked pixels

        return noise


