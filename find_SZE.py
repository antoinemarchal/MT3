from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

import function as fct
plt.ion()

patch_size = 100 #FIXME
PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT = fct.coord_sz(PSZ)
#w.wcs.crval = [121.1743 ,-21.5733]

new_map = [] #declare une liste

unit_1 = open("filenames_HFI.txt")

for line in unit_1:
    filename_smooth = line.strip()
    map_smooth,header = hp.read_map(filename_smooth,h=True)
    new_map.append(fct.patch_map(map_smooth, patch_size, GLON[0], GLAT[0]))
    
plt.imshow(new_map[0]) #voir dictionnaire



