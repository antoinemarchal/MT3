from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

import function as fct
plt.ion()

patch_size = 200 
PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT = fct.coord_SZ(PSZ)
#test for andromeda galaxy [121.1743 ,-21.5733]

freq          = [100., 143., 217., 353., 545., 857.]
a             = np.ones(patch_size) #FIXME
a_t           = np.transpose(a)     #FIXME
b             = []
b_t           = []
patch_map     = []

w_t           = []
f_nu          = fct.fact_SZE(freq)
print f_nu

unit_1 = open("files_HFI_full.txt")
path_1 = "maps_smooth/"
i = 0
for line in unit_1:
    filename_smooth = line.strip()
    map_smooth,header = hp.read_map(path_1 + filename_smooth[10:],h=True)
    patch_map.append(
        (fct.patch_map(map_smooth, patch_size, GLON[0], GLAT[0]))
        )

    b.append([f_nu[i]] * patch_size) 
    b_t.append(np.transpose(b[i]))  
    
    i += 1

#FIXME Covariance between all maps 
E = np.cov((patch_map[0].flatten(), patch_map[1].flatten(),
              patch_map[2].flatten(), patch_map[3].flatten(),
              patch_map[4].flatten(), patch_map[5].flatten()))
C = np.linalg.inv(E)
C = np.sum(C) #FIXME

for i in range(5):
    w_t.append( (((a_t*C*a)*b_t[i]*C)-     #FIXME nan
              ((b_t[i]*C*a)*a_t*C)) /      
              ((b_t[i]*C*b[i])*(a_t*C*a)-(b_t[i]*C*a)**2)
    )
    print w_t[i] 

CMB_KSZ_map   = np.zeros((patch_size, patch_size))
TSZ_map = np.zeros((patch_size, patch_size))

for i in range(5) :
#    CMB_KSZ_map = CMB_KSZ_map + (w_k[i] * patch_map[i])
    TSZ_map     = TSZ_map     + (w_t[i] * patch_map[i])

#plt.imshow(TSZ_map)




