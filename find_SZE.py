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

patch_size = 100 #FIXME
"""------------------------------------------------------------
   -------Reading galactic coordonate into .fits file ---------
   ------------------------------------------------------------"""
PSZ = "PSZ2v1.fits"
GLON,GLAT = fct.coord_sz(PSZ)
#w.wcs.crval = [121.1743 ,-21.5733]

new_map = fct.patch_map(map_smooth, patch_size, GLON, GLAT, 0)

plt.imshow(new_map)



