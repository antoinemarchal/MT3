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
plt.ion()

n_cluster      = 20
patch_size     = 200 
PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT = fct.coord_SZ(PSZ)

st_w          = []

n_id          = 20
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
    CMB_KSZ_map, TSZ_map, w_t, w_k = fct.separation(k, unit_1, origin
                                , patch_size, NAME, GLON, GLAT, freq)
    fct.save_fits(NAME[k], TSZ_map,k)
    sys.stdout = sys.__stdout__

    print k+1, " of ", n_cluster
    st_w.append(w_t[0])

plt.figure()
bins = np.linspace(-0.1, 0.1, 100)
plt.hist(st_w, bins)
plt.figure()
plt.plot(GLAT[:n_cluster],st_w, ".")

with open('output.dat', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(st_w)

sys.exit('stop')
#"""---------------Plot CMB+KSZ and TSZ--------------"""
plt.suptitle(NAME[22], size=16)

plt.subplot(1,2,1)
plt.imshow(TSZ_map, origin='lower', interpolation='none')
plt.title('TSZ map')

plt.subplot(1,2,2)
plt.imshow(CMB_KSZ_map, origin='lower', interpolation='none')
plt.title('CMB+KSZ map')
        
    



