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
moy_w = []
for i in range(n_obs):
    moy_w.append(np.mean(w[:n_cluster,i]))

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
"""---------Creating patches of SZ effect--------
-------------------------------------------------"""    
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
    RD = []
    for i in range(n_obs):
        RD.append(
            np.absolute((w_t[i]-moy_w[i])) / moy_w[i]
        )
    #if k == 22:
    #    inout.plot_map(NAME[k], TSZ_map, CMB_KSZ_map)
    
    if (np.absolute(RD[0]) > 0.90):
        GLAT_excl[j]   = GLAT[k]
        st_w_excl[j,:] = w_t
        j += 1
    else:
        st_w[l,:]    = w_t
        GLAT_slct[l] = GLAT[k]
        #inout.save_fits(NAME[k], TSZ_map,k)
        l += 1

    print k+1#, ':', k+1-j, 'of', n_cluster, 'selected'

n_slct = l
n_excl = j

#with open('sort_weight.dat', 'wb') as output:
#    mon_pickler = pickle.Pickler(output)
#    mon_pickler.dump(st_w)

inout.plot_weight(GLAT_slct, GLAT_excl, st_w, st_w_excl, n_slct, n_excl)
