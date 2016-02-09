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

n_cluster      = 30#1653
patch_size     = 100 
PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT = inout.coord_SZ(PSZ)

n_obs         = 6
st_w          = np.zeros((n_cluster,n_obs) )
fwhm          = np.loadtxt("fwhm_HFI.txt")
freq          = []
for i in range(n_obs):
    freq.append(fwhm[i,0])

"""---Compute the average value of weight for 
      each frequency-----------------------------"""
w = pickle.load(open("weight.dat","rb"))
moy_w = []
for i in range(6):
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
            np.absolute((w_t[i]-moy_w[i]) / moy_w[i])
        )
    if k == 22:
        inout.plot_map(NAME[k], TSZ_map, CMB_KSZ_map)
    
    #if (np.amin(RD) < 0.1):
    st_w[j,:] = w_t
    #inout.save_fits(NAME[k], TSZ_map,k)
    j += 1
        
    print k+1, ':', j, 'of', n_cluster, 'selected'
    
n_slct = j

with open('sort_weight.dat', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(st_w)

inout.plot_weight(GLAT, st_w, n_slct)
