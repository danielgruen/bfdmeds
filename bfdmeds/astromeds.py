# Derive a class from MEDS which overrides the astrometric solution with a new one
# found in PixelMapCollection YAML file.  We do NOT attempt to fix up the
# map for the cutout 0, which is always the coadd cutout.

import meds
from pixmappy import PixelMapCollection
import numpy as np
import re
import astropy.coordinates as co

# Map from CCD number to DETPOS value
detpos  =  {1:'S29' , 2:'S30', 3:'S31', 4 :'S25', 5 :'S26', 6 :'S27',
            7 :'S28', 8 :'S20', 9 :'S21', 10:'S22', 11:'S23', 12:'S24',
            13:'S14', 14:'S15', 15:'S16', 16:'S17', 17:'S18', 18:'S19',
            19:'S8' , 20:'S9' , 21:'S10', 22:'S11', 23:'S12', 24:'S13',
            25:'S1' , 26:'S2' , 27:'S3' , 28:'S4' , 29:'S5' , 30:'S6' ,
            31:'S7' , 32:'N1' , 33:'N2' , 34:'N3' , 35:'N4' , 36:'N5' ,
            37:'N6' , 38:'N7' , 39:'N8' , 40:'N9' , 41:'N10', 42:'N11',
            43:'N12', 44:'N13', 45:'N14', 46:'N15', 47:'N16', 48:'N17',
            49:'N18', 50:'N19', 51:'N20', 52:'N21', 53:'N22', 54:'N23',
            55:'N24', 56:'N25', 57:'N26', 58:'N27', 59:'N28', 60:'N29',
            61:'N30', 62:'N31'}

class AstroMEDS(meds.MEDS):
    def __init__(self, meds_file, astro_file=None):
        """
        Create MEDS file wrapped by new astrometry.

        parameters
        ----------
        meds_file:
            Filename for the MEDS catalog
        astro_file:
            Filename of YAML-format astrometric solutions.  None (default) to
            keep using astrometry in the MEDS file
        """
        super(AstroMEDS, self).__init__(meds_file)
        if astro_file is None:
            self.pmc = None
        else:
            self.pmc = PixelMapCollection(astro_file)
        self.color = 0.5 # ??? Nominal color given to all WCS's
        self.step = 5    # Step size (in pixels) used for Jacobian calcs

    def _get_wcs(self, file_id):
        '''Internal to retrieve the WCS for a given file
        '''
        file_re = re.compile(r'.*D(\d*)_._c(\d*).*')
        matches = file_re.match(file_id)
        if matches is None or len(matches.groups())!=2:
            print 'Error pulling expnum & ccd from file_id', file_id
            return None
        expnum = int(matches.group(1))
        ccdnum = int(matches.group(2))
        if ccdnum not in detpos:
            print 'Error interpreting ccdnum', ccdnum
            return None
        wcsname = 'D{:06d}/{:s}'.format(expnum,detpos[ccdnum])
        return self.pmc.getWCS(wcsname)

    def _old_orig_rowcol(self,iobj,icutout):
        # MEDS will eventually have a get_orig_rowcol() method, but doesn't
        # yet.  When it does this is simply
        # return super(AstroMEDS,self).get_orig_rowcol(iobj,icutout)
        return self['orig_row'][iobj][icutout], self['orig_col'][iobj][icutout]        

    def _old_cutout_rowcol(self,iobj,icutout):
        # Same thing for cutout rowcol
        return super(uberMEDS,self).get_cutout_rowcol(iobj,icutout)
        #return self['cutout_row'][iobj][icutout], self['cutout_col'][iobj][icutout]
    
    def get_orig_rowcol(self,iobj, icutout):
        """
        get orig_row, orig_col for the specified object
        and epoch, using PixelMapCollection.

        parameters
        ----------
        iobj:
            Index of the object
        icutout: integer
            Index of the cutout for this object.

        returns
        -------
        row,col the location in the cutout image
        """
        # Get old position as guess at new one
        old_rowcol = self._old_orig_rowcol(iobj,icutout)
        if icutout==0:
            # No change for coadd cutout:
            return old_rowcol
        ra,dec = self['ra'][iobj],self['dec'][iobj]
        posn = co.SkyCoord(ra, dec, unit='deg',frame='icrs')
        file_id = self.get_image_info()['image_path'][self['file_id'][iobj][icutout]]
        wcs = self._get_wcs(file_id)
        
        # Old value is guess for solver for new one. Swap rowcol order for xy
        guessxy = np.array( [old_rowcol[1], old_rowcol[0]])
        # Now solve.
        xy = wcs.toPix(posn, c=self.color, guess=guessxy)
        return xy[1],xy[0]   # rowcol is reversed order

    def get_cutout_rowcol(self,iobj, icutout):
        """
        get cutout_row, cutout_col for the specified object
        and epoch, using PixelMapCollection.

        parameters
        ----------
        iobj:
            Index of the object
        icutout: integer
            Index of the cutout for this object.

        returns
        -------
        row,col the location in the cutout image
        """
        if icutout==0 or self.pmc is None:
            # Use old version for coadd cutout
            return super(AstroMEDS,self).get_cutout_rowcol(iobj,icutout)
        # Subtract cutout corner from original coordinates to get cutout coords
        row,col = self.get_orig_rowcol(iobj, icutout)
        return row - self['orig_start_row'][iobj][icutout], \
               col - self['orig_start_col'][iobj][icutout]

    def get_coutout_rowcol_list(self,iobj):
        if self.pmc is None:
            return super(AstroMEDS,self).get_cutout_rowcol_list(iobj)
        ra,dec = self['ra'][iobj],self['dec'][iobj]
        posn = co.SkyCoord(ra, dec, unit='deg',frame='icrs')
        nstamps = self['ncutout'][iobj]
        file_id = self.get_image_info()['image_path'][self['file_id'][iobj][1:nstamps]]
        row0 = m['orig_start_row'][index][1:nstamps]
        col0 = m['orig_start_col'][index][1:nstamps]
        oldrow = m['cutout_row'][index][1:nstamps]
        oldcol = m['cutout_col'][index][1:nstamps]

        out = [ self._old_cutout_rowcol(iobj,0) ] # Pass coadd through unaltered
        for i,f in enumerate(file_id):
            guessxy=np.array((oldcol[i],oldrow[i]))
            wcs = self._get_wcs(f)
            xy = wcs.toPix(posn, c=self.color, guess=guessxy)
            out.append((xy[1] - row0[i],xy[0]-col0[i])) # "rowcol" is reverse of "xy"
        return out

    def get_jacobian(self,iobj, icutout):
        """
        Get the jacobian as a dict keyed by

            row0
            col0
            dudrow
            dudcol
            dvdcol
            dvdrow

        using PixelMapCollection.

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of the cutout for this object.
        """
        if icutout==0 or self.pmc is None:
            # No change for coadd cutout:
            return super(AstroMEDS,self).get_jacobian(iobj,icutout)
        
        # Get old position as guess at new one
        old_rowcol = self._old_orig_rowcol(iobj,icutout)
        ra,dec = self['ra'][iobj],self['dec'][iobj]
        posn = co.SkyCoord(ra, dec, unit='deg',frame='icrs')
        file_id = self.get_image_info()['image_path'][self['file_id'][iobj][icutout]]
        wcs = self._get_wcs(file_id)
        
        # Old value is guess for solver for new one. Swap rowcol order for xy
        guessxy = np.array( [old_rowcol[1], old_rowcol[0]])
        # Now solve.
        xy = wcs.toPix(posn, c=self.color, guess=guessxy)
        # And get Jacobian
        jac = wcs.jacobian(xy, c=self.color, step=self.step)
        # Convert sky units from degrees to arcsec
        jac *= 3600.
        
        result = {'row0':(xy[1] - self['orig_start_row'][iobj][icutout]),
                  'col0':(xy[0] - self['orig_start_col'][iobj][icutout]),
                  'dudrow':-jac[0][1],
                  'dvdrow':jac[1][1],
                  'dudcol':-jac[0][0],
                  'dvdcol':jac[1][0]}
        
        return result
        
    def get_jacobian_matrix(self,iobj, icutout):
        """
        Get the jacobian as a numpy matrix, using PixelMapCollection.


        parameters
        ----------
        iobj:
            Index of the object
        icutout: integer
            Index of the cutout for this object.

        returns
        -------
        A 2x2 matrix of the jacobian
            dudrow dudcol
            dvdrow dvdcol
        """
        if icutout==0 or self.pmc is None:
            # No change for coadd cutout:
            return super(AstroMEDS,self).get_jacobian_matrix(iobj,icutout)
        
        d = self.get_jacobian(iobj, icutout)
        return  np.array( [d['dudrow'], d['dudcol']],
                          [d['dvdrow'], d['dvdcol']] )
        
    def get_jacobian_list(self,iobj):
        if self.pmc is None:
            return super(AstroMEDS,self).get_jacobian_list(iobj)
        nstamps = self['ncutout'][iobj]
        return [self.get_jacobian(iobj,icutout) for icoutout in range(nstamps)]
    

        

