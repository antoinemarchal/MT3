from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

import function as fct
plt.ion()

patch_size = 400 
PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT = fct.coord_SZ(PSZ)
#test for andromeda galaxy [121.1743 ,-21.5733]

n_obs         = 6

freq          = [100., 143., 217., 353., 545., 857.]
a             = np.ones(n_obs)
a_t           = np.transpose(a)     
patch_map     = []

w_t           = []
f_nu          = fct.fact_SZE(freq)
b             = np.ones(n_obs)
b             = f_nu
#for i in range(n_obs):
#    b[i] = b[i] * f_nu[i]
b_t           = np.transpose(b)

unit_1 = open("files_HFI_full.txt")
path_1 = "maps_smooth/"
i = 0
for line in unit_1:
    filename_smooth = line.strip()
    map_smooth,header = hp.read_map(path_1 + filename_smooth[10:],h=True)
    patch_map.append(
        (fct.patch_map(map_smooth, patch_size, GLON[0], GLAT[0]))
        ) 
    i += 1

#FIXME Covariance between all maps 
E = np.cov((patch_map[0].flatten(), patch_map[1].flatten(),
              patch_map[2].flatten(), patch_map[3].flatten(),
              patch_map[4].flatten(), patch_map[5].flatten()))
R = np.linalg.inv(E)

f_1 = np.dot(np.dot(a_t,R),a) * (np.dot(b_t,R))
f_2 = np.dot(np.dot(b_t,R),a) * (np.dot(a_t,R))
f_3 = np.dot(np.dot(b_t,R),b) * np.dot(np.dot(a_t,R),a)
f_4 = np.dot(np.dot(b_t,R),a)**2

w_t = (f_1 - f_2) / (f_3 - f_4)
print w_t

#    w_t.append( (((a_t[i]*C*a[i])*b_t[i]*C)-     #FIXME nan
#              ((b_t[i]*C*a[i])*a_t[i]*C)) /      
#              ((b_t[i]*C*b[i])*(a_t[i]*C*a[i])-(b_t[i]*C*a[i])**2)
#    )
#w_t = ((np.dot(a_t,C,a)*np.dot(b_t,C))-((b_t*C*a)*a_t*C)) / \
#      ((b_t*C*b)*(a_t*C*a)-(b_t*C*a)**2)

CMB_KSZ_map   = np.zeros((patch_size, patch_size))
TSZ_map       = np.zeros((patch_size, patch_size))

for i in range(5) :
#    CMB_KSZ_map = CMB_KSZ_map + (w_k[i] * patch_map[i])
    TSZ_map     = TSZ_map     + (w_t[i] * patch_map[i])

plt.imshow(TSZ_map)




