from astropy.io import fits as pyfits
import pyfits as pf
from astropy import wcs
import astropy.table as pytabs
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os
import function as fct
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

def save_fits(name, patch,indice):
    path = "patch_SZ/SZ/"
    filename = path + str(indice)+'_' + name + ".fits"
    if os.path.isfile(filename) == 1:
        os.remove(filename)    
    pf.writeto(filename, patch.value)
    return 0

def plot_map(NAME, TSZ_map, CMB_KSZ_map):
    #"""---------------Plot CMB+KSZ and TSZ--------------"""
    plt.suptitle(NAME, size=16)

    plt.subplot(1,2,1)
    plt.imshow(TSZ_map, origin='lower', interpolation='none')
    plt.title('TSZ map')
    
    plt.subplot(1,2,2)
    plt.imshow(CMB_KSZ_map, origin='lower', interpolation='none')
    plt.title('CMB+KSZ map')
    return 0

def plot_weight(GLAT, st_w, n_cluster):
    #st_w = pickle.load(open("output.dat","rb"))
    plt.figure()
    bins = np.linspace(-0.1, 0.1, 100)
    plt.hist(st_w, bins)
    plt.figure()
    plt.plot(GLAT[:n_cluster],st_w, ".")
    return 0 
