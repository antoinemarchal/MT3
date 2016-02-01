from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

import function as fct
plt.ion()

patch_size = 400 #FIXME
PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT = fct.coord_SZ(PSZ)
#test for andromeda galaxy [121.1743 ,-21.5733]

a             = np.ones(patch_size)
a_t           = np.transpose(a)
SZ_map        = []
patch_map     = []
inv_cov       = []
freq          = [100, 143, 217, 353, 545, 857]

unit_1 = open("files_HFI_full.txt")
path_1 = "maps_smooth/"
i = 0
for line in unit_1:
    filename_smooth = line.strip()
    map_smooth,header = hp.read_map(path_1 + filename_smooth[10:],h=True)
    patch_map.append(
        (fct.patch_map(map_smooth, patch_size, GLON[0], GLAT[0]))
        )
    inv_cov.append(np.linalg.inv(np.cov(patch_map[i])))
    SZ_map.append((a_t * inv_cov[i] / a_t / inv_cov[i] / a) * patch_map[i])
    i += 1
    
#plt.imshow(SZ_map[4])
#plt.imshow(inv_cov[0])

#f_nu = fct.fact_SZE(freq) 

final_map = np.zeros((patch_size, patch_size))
for i in range(5) :
    final_map = final_map + SZ_map[i]

#plt.imshow(final_map)
