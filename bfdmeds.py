"""
Defines a derived class of MEDS (Multi Epoch Data Structures) to provide the
input BFD needs.

TODO: make sure Katie and I agree on an interface
TODO: implement the getters properly
TODO: figure out at which point we determine which exposures to cull (because of masking / non-white noise) and how that culling happens

"""

import meds
import numpy as np
import scipy.stats

class BFDMEDS(meds.MEDS):

    def get_cutout_list(self, iobj, type='image'):
        """
        Get an image list with all cutouts associated with this coadd object,
        subtracting neighbors and filling missing pixels in each cutout with 
        the MOF models 

        Note each individual cutout is actually a view into a larger
        mosaic of all images.

        parameters
        ----------
        iobj:
            Index of the object
        type: string, optional
            Cutout type. Default is 'image'.  Allowed
            values are 'image','weight','seg' 
            (the latter two are taken from the base class as is)

        returns
        -------
        A list of images hold all cutouts.
        """

        if(type!='image'):
          return super(BFDMEDS, self).get_cutout_list(iobj, type)

        # TODO: fill bad pixels instead
        return super(BFDMEDS, self).get_cutout_list(iobj, 'image')


    def get_rowcol_list(self, iobj):
        """
        Get a list of centroid rows and columns for all cutouts associated with this
        coadd object.

        TODO: decide whether we need this since the same information is contained in
              the output of get_jacobian_list

        parameters
        ----------
        iobj:
            Index of the object

        returns
        -------
        A list of row positions of object centroid, and a list of its col positions.
        """

        rowlist=[]
        collist=[]

        ncutout=self._cat['ncutout'][iobj]

        for i in xrange(ncutout):
            row,col = self.get_cutout_rowcol(iobj,i)
            rowlist.append(row)
            collist.append(col)

        return rowlist, collist


    def get_psf_list(self, iobj):
        """
        Get a list of PSF postage stamp images for all cutouts associated with this
        coadd object.

        parameters
        ----------
        iobj:
            Index of the object

        returns
        -------
        A list of PSF postage stamps.
        """

        ncutout=self._cat['ncutout'][iobj]

        psflist=[]
        for i in xrange(ncutout):
            psflist.append("") # TODO: implement

     
    # def get_jacobian_list(self) is already implemented with all we need in the base class


    def get_noise_list(self, iobj):
        """
        Get a list of mean rms noise levels for all cutouts associated with this
        coadd object, which is the inverse square root of the harmonic mean of weights.
        """

        ncutout=self._cat['ncutout'][iobj]

        wimgs=self.get_cutout_list(iobj,type='weight')

        noise=[]

        for i in xrange(ncutout):
            noise.append(1./np.sqrt(scipy.stats.hmean(wimgs[i][wimgs[i]>0],axis=None)))

        # TODO: deal with bad/masked pixels

        return noise


