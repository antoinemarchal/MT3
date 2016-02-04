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

n_id          = 22
n_obs         = 6
FWHM          = np.loadtxt("FWMH_HFI.txt")
freq          = []
for i in range(n_obs):
    freq.append(FWHM[i,0])

patch_map     = []

w_t           = []
f_nu          = fct.fact_SZE(freq)

a             = np.ones(n_obs)
a_t           = np.transpose(a)
b             = np.ones(n_obs)
b             = f_nu
b_t           = np.transpose(b)

unit_1 = open("files_HFI_full.txt")
path_1 = "maps_smooth/"
i = 0
for line in unit_1:
    filename_smooth = line.strip()
    map_smooth,header = hp.read_map(path_1 + filename_smooth[10:],h=True)
    patch_map.append(
        (fct.patch_map(map_smooth, patch_size, GLON[n_id], GLAT[n_id]))
        ) 
    i += 1

#FIXME with n_obs
E = np.cov((patch_map[0].flatten(), patch_map[1].flatten(),
              patch_map[2].flatten(), patch_map[3].flatten(),
              patch_map[4].flatten(), patch_map[5].flatten()))
R = np.linalg.inv(E)

(w_k, w_t) = fct.weight(R, a, a_t, b ,b_t)

CMB_KSZ_map   = np.zeros((patch_size, patch_size))
TSZ_map       = np.zeros((patch_size, patch_size))

for i in range(n_obs) :
    CMB_KSZ_map = CMB_KSZ_map + (w_k[i] * patch_map[i])
    TSZ_map     = TSZ_map     + (w_t[i] * patch_map[i])

"""---------------------Plot CMB+KSZ and TSZ-------------------------"""
plt.suptitle(NAME[n_id], size=16)

plt.subplot(1,2,1)
plt.imshow(TSZ_map)
plt.title('TSZ map')

plt.subplot(1,2,2)
plt.imshow(CMB_KSZ_map)
plt.title('CMB+KSZ map')



