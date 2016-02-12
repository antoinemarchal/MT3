from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import pickle
import sys
import os

import mod_ilc as ilc
import in_output as inout

plt.ion()

n_cluster      = 1653
patch_size     = 100 
PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT = inout.coord_SZ(PSZ)

n_obs         = 6
st_w          = np.zeros((n_cluster,n_obs))
GLAT_slct     = np.zeros(n_cluster)
st_w_excl     = np.zeros((n_cluster,n_obs))
GLAT_excl     = np.zeros(n_cluster)
fwhm          = np.loadtxt("fwhm_HFI.txt")
freq          = []
for i in range(n_obs):
    freq.append(fwhm[i,0])

"""---Compute the average value of weight for 
      each frequency-----------------------------"""
w = pickle.load(open("weight.dat","rb"))
moy_w      = []
moy_w_slct = []
std_w_slct = []
std_w      = []
for i in range(n_obs):
    moy_w.append(np.mean(w[:n_cluster,i]))
    std_w.append(np.std(w[:n_cluster,i]))

print "********************************************"
print "********************************************"
print "********************************************"

"""------------Reading the smooth maps -------------
-------------------------------------------------"""
unit_1 = open("files_HFI_full.txt")
#unit_1.seek(origin)
#origin = unit_1.tell()
path_1 = "maps_smooth/"
map_smooth = []
header     = []
for line in unit_1:
    sys.stdout = open(os.devnull, "w")
    filename_smooth = line.strip()
    map_smooth.append(
        hp.read_map(path_1 + filename_smooth)
    )
"""-----------Separation of full sky map------------
-------------------------------------------------"""
f_nu          = ilc.dist_SZ(freq)
a             = np.ones(n_obs)
a_t           = np.transpose(a)
b             = np.ones(n_obs)
b             = f_nu
b_t           = np.transpose(b)
CMB   = np.zeros(hp.get_map_size(map_smooth[0]))
TSZ   = np.zeros(hp.get_map_size(map_smooth[0]))
J = np.cov((map_smooth[0], map_smooth[1], map_smooth[2], map_smooth[3],
           map_smooth[4], map_smooth[5]))
K = np.linalg.inv(J)
WK, WT = ilc.weight(K, a, a_t, b ,b_t)
for i in range(n_obs) :
    CMB = CMB + (WK[i] * map_smooth[i])
    TSZ = TSZ + (WT[i] * map_smooth[i])

"""---------Creating patches of SZ effect--------
   ----------------------------------------------"""    
print "Creating patches of SZ effect : "
j = 0
l = 0
for k in range(n_cluster):
    sys.stdout = open(os.devnull, "w")
    CMB_KSZ_map, TSZ_map, w_t, w_k = ilc.separation(
        k, map_smooth, patch_size, NAME, GLON,
        GLAT, freq
    )
    sys.stdout = sys.__stdout__
    
    """------Separation of cluster population------
       --------------------------------------------"""
    RD = []
    for i in range(n_obs-2):
        RD.append(
            np.absolute((w_t[i]-moy_w[i]) / std_w[i])
        )
    if k == 22:
        inout.plot_map(NAME[k], TSZ_map, CMB_KSZ_map)

    if ((RD[0] > 1.5)
    | (RD[1] > 1.5)
    | (RD[2] > 1.5)):
        GLAT_excl[j]   = GLAT[k]
        st_w_excl[j,:] = w_t
        j += 1
    else:
        st_w[l,:]      = w_t
        GLAT_slct[l]   = GLAT[k]
        #inout.save_fits(NAME[k], TSZ_map,k)
        l += 1
        
    print k+1#, ':', k+1-j, 'of', n_cluster, 'selected'
n_slct = l
n_excl = j

for i in range(n_obs):
    moy_w_slct.append(np.mean(st_w[:n_slct,i]))
    std_w_slct.append(np.std(st_w[:n_slct,i]))

w_full_sky = np.loadtxt("w_full_sky.d")

inout.plot_w_full_patch(w_full_sky[:,0], moy_w_slct, moy_w, std_w_slct,
                        std_w, WT)
inout.plot_weight(GLAT_slct, GLAT_excl, st_w, st_w_excl, n_slct, n_excl)



    
