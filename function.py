from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

plt.ion()

def smooth(filename, fwhm) :
        map,header = hp.read_map(filename,h=True)
        # map_size = hp.get_map_size(map)
        #hp.mollview(map, norm='hist')
        # print(hp.npix2nside(map_size))
        
        print fwhm
        fwhm_rad = np.radians(fwhm/60.)
        print fwhm_rad
        map_gauss = hp.smoothing(map, fwhm_rad)
        #hp.mollview(map_gauss, norm='hist')
        print 'end of smooth'

        return (map_gauss, header)

def coord_sz(filename) :
        import astropy.table as pytabs
        #import astropy.io.fits as af
        
        cat  = pyfits.getdata(filename)
        data = pytabs.Table(cat)
        NAME = data['NAME']
        RA   = data['RA']
        DEC  = data['DEC']

        return (NAME,RA,DEC)
    
def patch_map(map_smooth, patch_size, GLON, GLAT) :
        """------------------------------------------------------------
        -------------Creation of wcs objects for patch--------------
        ------------------------------------------------------------"""
        w           = wcs.WCS(naxis=2)
        w.wcs.crpix = [patch_size/2, patch_size/2] #FIXME +1 si 0 ou 1
        w.wcs.crval = [GLON,GLAT]
        n_side      = hp.get_nside(map_smooth)
        pix_size    = ma.sqrt(hp.nside2pixarea(n_side, degrees=True))
        w.wcs.cdelt = np.array([-pix_size/2., pix_size/2.])
        w.wcs.ctype = ["GLON-TAN", "GLAT-TAN"]     # use a gnomonic projetion
    
        header = w.to_header() #FIXME
    
        patch   = np.zeros((patch_size, patch_size))
        new_map = patch
        xx, yy  = np.indices((patch_size, patch_size))
    
        patch_GLON, patch_GLAT = w.wcs_pix2world(xx, yy, 0)
        patch_PHI = patch_GLON * ma.pi / 180. 
        patch_THETA = patch_GLAT * ma.pi / 180.
        patch = hp.ang2pix(n_side,ma.pi / 2. - patch_THETA, patch_PHI)
        new_map = map_smooth[patch]

        return new_map
