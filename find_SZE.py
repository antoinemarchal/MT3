from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

import function as fct
plt.ion()

patch_size = 2000 #FIXME
PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT = fct.coord_SZ(PSZ)
#test for andromeda galaxy [121.1743 ,-21.5733]

patch_map = [] 
#patch_map = {}

unit_1 = open("filenames_HFI.txt")
path_1 = "maps_smooth/"
for line in unit_1:
    filename_smooth = line.strip()
    map_smooth,header = hp.read_map(path_1 + filename_smooth[10:],h=True)
   patch_map.append(
        (filename_smooth[21:24], fct.patch_map(map_smooth, patch_size, GLON[0], GLAT[0]))
        )
#plt.imshow(patch_map[0][1]) 
#    patch_map["" + filename_smooth[21:24]] = fct.patch_map(map_smooth, patch_size, GLON[0], GLAT[0])
#plt.imshow(patch_map[0])
