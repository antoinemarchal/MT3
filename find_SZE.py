from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

import fonction as fct

plt.ion()
filename_smooth = "maps_smooth/HFI_SkyMap_857_2048_R2.00_full.fits"
map_smooth,header = hp.read_map(filename_smooth,h=True)
#hp.mollview(map_smooth, norm='hist')
"""
Creation of wcs objects for patch

FIXME information

"""
PSZ = "PSZ2v1.fits"
GLON,GLAT = fct.coord_sz(PSZ)

patch_size = 200 #FIXME

w = wcs.WCS(naxis=2)
w.wcs.crpix = [patch_size/2, patch_size/2] #FIXME +1 si 0 ou 1
#w.wcs.crval = [GLON[0],GLAT[0]]
w.wcs.crval = [121.1743 ,-21.5733]
#map_size = hp.get_map_size(map_smooth)
n_side   = hp.get_nside(map_smooth)
#pix_size = ma.sqrt(4 * ma.pi / map_size) * 180. / ma.pi
pix_size = ma.sqrt(hp.nside2pixarea(n_side, degrees=True))
w.wcs.cdelt = np.array([-pix_size/2., pix_size/2.])        #FIXME
# use a gnomonic projetion
w.wcs.ctype = ["GLON-TAN", "GLAT-TAN"]

header = w.to_header()

patch   = np.zeros((patch_size, patch_size))
new_map = patch
xx, yy  = np.indices((patch_size, patch_size))

patch_GLON, patch_GLAT = w.wcs_pix2world(xx, yy, 0)
patch_PHI = patch_GLON * ma.pi / 180. 
patch_THETA = patch_GLAT * ma.pi / 180.
patch = hp.ang2pix(n_side,ma.pi / 2. - patch_THETA, patch_PHI)
print patch
new_map = map_smooth[patch]
print new_map
plt.imshow(new_map)



