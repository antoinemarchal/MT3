from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

import fonction as fct

plt.ion()
filename_smooth = "maps_smooth/HFI_SkyMap_143_2048_R2.00_full.fits"
map_smooth,header = hp.read_map(filename_smooth,h=True)
#hp.mollview(map_smooth, norm='hist')
"""
Creation of wcs objects for patch

FIXME information

"""
PSZ = "PSZ2v1.fits"
RA,DEC = fct.coord_sz(PSZ)

p_size = 200 #FIXME

w = wcs.WCS(naxis=2)
w.wcs.crpix = [p_size/2, p_size/2] #FIXME +1 si 0 ou 1
w.wcs.crval = [RA[0],DEC[0]]
#map_size = hp.get_map_size(map_smooth)
n_side   = hp.get_nside(map_smooth)
#pix_size = ma.sqrt(4 * ma.pi / map_size) * 180. / ma.pi
pix_size = ma.sqrt(hp.nside2pixarea(n_side, degrees=True))
w.wcs.cdelt = np.array([-pix_size/2., pix_size/2.])        #FIXME
# use a gnomonic projetion
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

header = w.to_header()

patch = np.zeros([p_size, p_size])
grid  = np.indices([p_size, p_size])
print grid
print grid.shape

p_RA, p_DEC = w.wcs_pix2world()

