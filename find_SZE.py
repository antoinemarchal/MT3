from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import pickle
import sys
import os

import function as fct
import mod_ILC as ilc
plt.ion()

n_cluster      = 2#1653
patch_size     = 200 
PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT = fct.coord_SZ(PSZ)

st_w          = []

n_obs         = 6
FWHM          = np.loadtxt("FWMH_HFI.txt")
freq          = []
for i in range(n_obs):
    freq.append(FWHM[i,0])

unit_1 = open("files_HFI_full.txt")
origin = unit_1.tell()

print "***********************************************************"
print "***********************************************************"
print "***********************************************************"

print "Creating patches of SZ effect : "
for k in range(n_cluster):
    sys.stdout = open(os.devnull, "w")
    CMB_KSZ_map, TSZ_map, w_t, w_k = ilc.separation(k, unit_1, origin
                                , patch_size, NAME, GLON, GLAT, freq)
    fct.save_fits(NAME[k], TSZ_map,k)
    sys.stdout = sys.__stdout__

    print k+1, " of ", n_cluster
    st_w.append(w_t[0])
    if k == 22:
        fct.plot_map(NAME[k], TSZ_map, CMB_KSZ_map)

#with open('weight_0.dat', 'wb') as output:
#    mon_pickler = pickle.Pickler(output)
#    mon_pickler.dump(st_w)

