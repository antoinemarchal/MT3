from astropy.io import fits as pyfits
from astropy import wcs
import astropy.table as pytabs
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os
plt.ion()

def smooth(filename, fwhm) :
        """-------------------------------------------------
        --- smooth : Convuluted form by a gaussian.
                        
                     Parameters : filename  = adress of map
                                     
                     Return     : map_gauss = convoluted map
                                  header
        ----------------------------------------------------"""
        map,header = hp.read_map(filename,h=True)
        fwhm_rad = np.radians(fwhm/60.)
        map_gauss = hp.smoothing(map, fwhm_rad)
        return (map_gauss, header)

def coord_SZ(filename) :
        """-------------------------------------------------
        --- coord_SZ : Getting back useful information of
                       galaxy cluster in SZ catalog 
                        
                       Parameters : filename  = .fits file
                                     
                       Return     : NAME = convoluted map
                                  GLON = galactic longitude
                                  GLAT = galactic latitude
        ----------------------------------------------------"""
        cat  = pyfits.getdata(filename)
        data = pytabs.Table(cat)
        NAME = data['NAME']
        GLON = data['GLON']
        GLAT = data['GLAT']
        return (NAME,GLON,GLAT)
    
def patch_map(map_smooth, patch_size, GLON, GLAT) :
        """------------------------------------------------------------
        --- patch_map : Initiate the WCS object and return a piece of
                        Planck map centred on a given coordonates
                        using a gnomonic projection.
                        
                        Parameters : GLON       = galactic longitude
                                     GLAT       = galactic latitude
                                     patch_size = #FIXME dimension

                        Return     : new_map (2D), D = patch_size
        ------------------------------------------------------------"""
        w           = wcs.WCS(naxis=2)
        w.wcs.crpix = [patch_size/2, patch_size/2] #FIXME +1 si 0 ou 1
        w.wcs.crval = [GLON,GLAT]
        n_side      = hp.get_nside(map_smooth)
        pix_size    = ma.sqrt(hp.nside2pixarea(n_side, degrees=True))
        w.wcs.cdelt = np.array([-pix_size/2., pix_size/2.])
        w.wcs.ctype = ["GLON-TAN", "GLAT-TAN"] 
    
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

    def fact_SZE(nu) :
        import astropy.units as u
        from astropy.constants import G, h, k_B
        from astropy.cosmology import FLRW
        x = u.h * nu / u.k_B / FLRW.Tcmb0
        f = x * ((ma.exp(x) + 1.) / (ma.exp(x) - 1.)) - 4.
        return f
