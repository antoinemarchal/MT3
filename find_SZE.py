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
freq      = [100, 143, 217, 353, 545, 857]

unit_1 = open("files_HFI_full.txt")
path_1 = "maps_smooth/"
for line in unit_1:
    filename_smooth = line.strip()
    map_smooth,header = hp.read_map(path_1 + filename_smooth[10:],h=True)
    patch_map.append(
        (fct.patch_map(map_smooth, patch_size, GLON[0], GLAT[0]))
        )
plt.imshow(patch_map[0]) 

f_nu = fct.fact_SZE(freq) 

    
