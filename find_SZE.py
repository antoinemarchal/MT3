from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

import function as fct
plt.ion()

filename_smooth = "maps_smooth/HFI_SkyMap_857_2048_R2.00_full.fits"
map_smooth,header = hp.read_map(filename_smooth,h=True)
#hp.mollview(map_smooth, norm='hist')

patch_size = 2000 #FIXME
"""------------------------------------------------------------
   -------Reading galactic coordonate into .fits file ---------
   ------------------------------------------------------------"""
PSZ = "PSZ2v1.fits"
GLON,GLAT = fct.coord_sz(PSZ)

"""------------------------------------------------------------
   -------------Creation of wcs objects for patch--------------
   ------------------------------------------------------------"""
w = wcs.WCS(naxis=2)
w.wcs.crpix = [patch_size/2, patch_size/2] #FIXME +1 si 0 ou 1
#w.wcs.crval = [GLON[0],GLAT[0]]
w.wcs.crval = [121.1743 ,-21.5733]
n_side      = hp.get_nside(map_smooth)
pix_size    = ma.sqrt(hp.nside2pixarea(n_side, degrees=True))
w.wcs.cdelt = np.array([-pix_size/2., pix_size/2.])
w.wcs.ctype = ["GLON-TAN", "GLAT-TAN"]               # use a gnomonic projetion

header = w.to_header() #FIXME

patch   = np.zeros((patch_size, patch_size))
new_map = patch
xx, yy  = np.indices((patch_size, patch_size))

patch_GLON, patch_GLAT = w.wcs_pix2world(xx, yy, 0)
patch_PHI = patch_GLON * ma.pi / 180. 
patch_THETA = patch_GLAT * ma.pi / 180.
patch = hp.ang2pix(n_side,ma.pi / 2. - patch_THETA, patch_PHI)
new_map = map_smooth[patch]

plt.imshow(new_map)



