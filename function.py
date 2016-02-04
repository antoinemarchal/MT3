from astropy.io import fits as pyfits
from astropy import wcs
import astropy.table as pytabs
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os
plt.ion()

def smooth(filename, fwhm):
        """-------------------------------------------------
        --- smooth : Convoluted form by a gaussian.
                        
                     Parameters : filename  = adress of map
                                     
                     Return     : map_gauss = convoluted map
                                  header
        ----------------------------------------------------"""
        map,header = hp.read_map(filename,h=True)
        fwhm_rad = np.radians(fwhm/60.)
        map_gauss = hp.smoothing(map, fwhm_rad)
        return (map_gauss, header)

def coord_SZ(filename):
        """--------------------------------------------------
        --- coord_SZ : Getting back useful information of
                       galaxy cluster in SZ catalog 
                        
                       Parameters : filename  = .fits file
                                     
                       Return     : NAME = ID of cluster
                                    GLON = galactic longitude
                                    GLAT = galactic latitude
        -----------------------------------------------------"""
        cat  = pyfits.getdata(filename)
        data = pytabs.Table(cat)
        NAME = data['NAME']
        GLON = data['GLON']
        GLAT = data['GLAT']
        return (NAME,GLON,GLAT)
    
def patch_map(map_smooth, patch_size, GLON, GLAT):
        """---------------------------------------------------------
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
        # Fortran or C convention ?
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
        patch_PHI              = patch_GLON * ma.pi / 180. 
        patch_THETA            = patch_GLAT * ma.pi / 180.
        patch                  = hp.ang2pix(n_side,ma.pi / 2. -
                                            patch_THETA, patch_PHI)
        new_map                = map_smooth[patch]
        return new_map

def fact_SZE(nu):
        """----------------------------------------------------
        --- fact_SZE : Getting back the factor f(nu) 
                        
                       Parameters : nu  = frequency in GHz
                                     
                       Return     : f   = adimenssionless factor
        --------------------------------------------------------"""
        import astropy.units as u
        from astropy import constants as const
        x = const.h * nu * u.GHz / const.k_B / ( 2.7 * u.K) #FIXME Tcmb
        #more precision
        f = x * ((np.exp(x) + 1.) / (np.exp(x) - 1.)) - 4.
        return f

def weight(R, a, a_t, b, b_t):
        """-------------------------------------------------------
        --- weight : Compute the weight of each map to obtein
                     CMB+KSZ map and TSZ map.
                        
                     Parameters : R      = matrix of covariance
                                           between the maps
                                  a, a_t = (1,...,1) and T matrix
                                  b, b_t = f(nu) and T matrix

                     Return     : w_k    = weight vector of CMB+KSZ 
                                           map
                                  w_t    = weight vector of TSZ map
        -----------------------------------------------------------"""
        g_1 = np.dot(np.dot(b_t,R),b) * (np.dot(a_t,R))
        g_2 = np.dot(np.dot(a_t,R),b) * (np.dot(b_t,R))
        g_3 = np.dot(np.dot(a_t,R),b) * np.dot(np.dot(b_t,R),b)
        g_4 = np.dot(np.dot(a_t,R),b)**2

        f_1 = np.dot(np.dot(a_t,R),a) * (np.dot(b_t,R))
        f_2 = np.dot(np.dot(b_t,R),a) * (np.dot(a_t,R))
        f_3 = np.dot(np.dot(b_t,R),b) * np.dot(np.dot(a_t,R),a)
        f_4 = np.dot(np.dot(b_t,R),a)**2

        w_k = (g_1 - g_2) / (g_3 - g_4)
        w_t = (f_1 - f_2) / (f_3 - f_4)
        return (w_k, w_t)


