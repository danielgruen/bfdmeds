"""
Defines a derived class of MEDS (Multi Epoch Data Structures) to provide the
input BFD needs.

TODO: make sure Katie and I agree on an interface
TODO: implement the getters properly
TODO: figure out at which point we determine which exposures to cull (because of masking / non-white noise) and how that culling happens

"""
import pdb
import meds
import numpy as np
import sys
import collections
from . import AstroMEDS

class BFDMEDS(AstroMEDS):


    def __init__(self, filename, psf_source=None, astro_dir=None, color_array=None, filename_mofsub=None, filename_mofsub_badfill=None):
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
            An instance of bfdpsf.PsfSource, if any. Default (None) says to use psf supplied in MEDS file
        astro_file:
            Filename of serialized WCSFit astrometry, if any.  Default (None) is to use MEDS data.
        filename_mofsub:
            MEDS filename, MOF neighbor models subtracted from all single epoch stamps
        fiename_mofsub_badfill:
            MEDS filename, MOF neighbor models subtracted and bad pixels filled
        """

        super(BFDMEDS, self).__init__(filename, astro_dir=astro_dir, color_array=color_array)        
        
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
            rows,cols = self.get_rowcol_list(iobj, False)
            return rows[1:], cols[1:]

        rowlist=[]
        collist=[]

        ncutout=self._cat['ncutout'][iobj]

        for i in xrange(ncutout):
            row,col = self.get_cutout_rowcol(iobj,i)
            rowlist.append(row)
            collist.append(col)

        return rowlist, collist


    def get_psf_list(self, iobj, skip_coadd=False,return_psf_object=False):
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
        if self.psf_source is not None:
            ncutout=self._cat['ncutout'][iobj]
            # make elegant to get coadd if needed
            
            if skip_coadd:
                return [self.get_psf(iobj, i, return_psf_object=return_psf_object) for i in xrange(1,ncutout)]
            else:
                return [self.get_psf(iobj, i, return_psf_object=return_psf_object) for i in xrange(0,ncutout)]
        else:
            ncut=self['ncutout'][iobj]
            psf_list=[super(BFDMEDS,self).get_psf(iobj,icut) for icut in xrange(ncut)]
            if skip_coadd:
                return psf_list[1:]
            else:
                return psf_list


    def get_psf(self, iobj, icutout,return_psf_object=False):
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
        
        row = self['orig_row'][iobj][icutout]
        col = self['orig_col'][iobj][icutout]

        stamp_size = self.get_cat()['box_size'][iobj]
        jacobian = self.get_jacobian(iobj,icutout)

        if icutout == 0:
            info = self.get_coadd_exposure_info(iobj)
            return self.psf_source.get_coadd_psf(info['tilename'], info['band'], info['ccd'], info['request_attempt'],col, row, stamp_size, jacobian,return_psf_object=return_psf_object)
 
        else:
            info = self.get_exposure_info(iobj, icutout)
            return self.psf_source.get_psf(info['tilename'], info['band'], info['exposure'], info['ccd'], col, row, stamp_size, jacobian,return_psf_object=return_psf_object)
            



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

        # trying to code this more foolproof
        tilename_start_index=image_path.find("DES")
        tilename=image_path[tilename_start_index:tilename_start_index+12]
        exposure_start_index=image_path.find("D00")
        exposure=image_path[exposure_start_index:exposure_start_index+9]
        band=image_path[exposure_start_index+10]
        ccd_part=image_path[exposure_start_index+12:exposure_start_index+15]
        request_attempt=image_path[exposure_start_index+16:exposure_start_index+24]

        #filename=path[-1]
        #tilename=path[12]        
        #exposure, band, ccd_part, request_attempt, _, _ = filename.split("_")

        #Strip out the boilerplate
        ccd = ccd_part.lstrip("c")

        #Return info as dictionary
        info = dict(ccd=ccd, tilename=tilename, request_attempt=request_attempt, band=band, exposure=exposure)
        return info

    def get_coadd_exposure_info(self, iobj):
        """
        Get a selection of information about this exposure.

        parameters
        ----------
        iobj:
            Index of the object

        returns
        -------
        A numpy array containing the PSF image.
        """
        # fix this for supernova fields vs. regular wide fields
        #Get the source info and thence image path
        info = self.get_source_info(iobj, 0)
        image_path = info['image_path'] 

        tile_start_index=image_path.find("DES")
        tilename=image_path[tile_start_index:tile_start_index+12]
        
        # image_paths have this format:
        # coadd/SN-C3_C28_r3499p02_r.fits.fz 
        # parse into encoded pieces
        coadd,filename=image_path.split("/")
        filename,fits=filename.split(".")
        tile,ccd,request_attempt,band=filename.split("_")

        #Return info as dictionary
        info = dict(ccd=ccd, tilename=tilename, request_attempt=request_attempt, band=band)
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
            return self.get_noise_list(iobj,skip_coadd=False)[1:]

        ncutout=self._cat['ncutout'][iobj]

        wimgs=self.get_cutout_list(iobj,type='weight')

        noise=[]

        for i in xrange(ncutout):
            noise.append(1./np.sqrt(np.median(wimgs[i][wimgs[i]>0],axis=None)))

        # TODO: deal with bad/masked pixels

        return noise


