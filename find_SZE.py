from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

import function as fct
plt.ion()

patch_size = 200 #FIXME
PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT = fct.coord_SZ(PSZ)
#test for andromeda galaxy [121.1743 ,-21.5733]

a             = np.ones(patch_size)
a_t           = np.transpose(a)

TSZ_map       = []
CMB_KSZ_map   = []
patch_map     = []
inv_cov       = []
freq          = [100, 143, 217, 353, 545, 857]

w             = []

f_nu          = fct.fact_SZE(freq)

unit_1 = open("files_HFI_full.txt")
path_1 = "maps_smooth/"
i = 0
for line in unit_1:
    filename_smooth = line.strip()
    map_smooth,header = hp.read_map(path_1 + filename_smooth[10:],h=True)
    patch_map.append(
        (fct.patch_map(map_smooth, patch_size, GLON[0], GLAT[0]))
        )

    b   = [f_nu[i]] * patch_size
    b_t = np.transpose(b)
    
    inv_cov.append(np.linalg.inv(np.cov(patch_map[i]))) #FIXME
    w.append( (((a_t*inv_cov[i]*a)*b_t*inv_cov[i])- 
              ((b_t*inv_cov[i]*a)*a_t*inv_cov[i])) / 
              ((b_t*inv_cov[i]*b)*(a_t*inv_cov[i]*a)-(b_t*inv_cov[i]*a)**2)
    )
    print w[i]
    
    #CMB_KSZ_map.append((a_t * inv_cov[i] / a_t / inv_cov[i] / a) \
    #* patch_map[i])
    
    #TSZ_map.append((b_t * inv_cov[i] / b_t / inv_cov[i] / b) \
    #* patch_map[i])
    i += 1
    
#plt.imshow(TSZ_map[0])
#plt.imshow(inv_cov[0])

#f_nu = fct.fact_SZE(freq) 

final_CMB_KSZ = np.zeros((patch_size, patch_size))
final_TSZ     = np.zeros((patch_size, patch_size))

for i in range(5) :
#    final_CMB_KSZ = final_CMB_KSZ + CMB_KSZ_map[i]
    final_TSZ     = final_TSZ     +  (w[i] * patch_map[i])

#for i in range(5) :
#    final_CMB_KSZ = final_CMB_KSZ + CMB_KSZ_map[i]
#    final_TSZ     = final_TSZ     + TSZ_map[i]

print final_TSZ
plt.imshow(final_TSZ)




