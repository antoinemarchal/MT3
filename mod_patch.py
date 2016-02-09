from astropy.io import fits as pyfits
import pyfits as pf
from astropy import wcs
import astropy.table as pytabs
import numpy as np
import healpy as hp
import math as ma

def patch_map(map_smooth, patch_size, GLON, GLAT):
    """---------------------------------------------------------
    --- patch_map : Initiate the WCS object and return a piece of
                    Planck map centred on a given coordonates
                    using a gnomonic projection.
    
                    Parameters : GLON       = galactic longitude
                                 GLAT       = galactic latitude 
                                 patch_size = number of pixel
    
                    Return     : new_map (2D), D = patch_size
    ------------------------------------------------------------"""
    w           = wcs.WCS(naxis=2)
    w.wcs.crpix = [patch_size/2, patch_size/2]
    w.wcs.crval = [GLON,GLAT]
    n_side      = hp.get_nside(map_smooth)
    pix_size    = ma.sqrt(hp.nside2pixarea(n_side, degrees=True))
    w.wcs.cdelt = np.array([-pix_size/2., pix_size/2.])
    w.wcs.ctype = ["GLON-TAN", "GLAT-TAN"] 
    
    header = w.to_header() #FIXME
    
    patch   = np.zeros((patch_size, patch_size))
    new_map = patch
    yy, xx  = np.indices((patch_size, patch_size))
    
    patch_GLON, patch_GLAT = w.wcs_pix2world(xx, yy, 0)
    patch_PHI              = patch_GLON * ma.pi / 180. 
    patch_THETA            = patch_GLAT * ma.pi / 180.
    patch                  = hp.ang2pix(n_side,ma.pi / 2. -
                                        patch_THETA, patch_PHI)
    new_map                = map_smooth[patch]
    return new_map

