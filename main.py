from astropy.io import fits as pyfits
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma

import fonction as fct
plt.ion()

unit_1 = open("filenames_HFI.txt")
FWMH = np.loadtxt("FWMH_HFI.txt")
i = 0
for line in unit_1 :
    fichier = "maps_2015/" + line.strip()
    if "100" in fichier :
        continue
    if "143" in fichier :
        j = 1
    if "217" in fichier :
        j = 2
    if "353" in fichier :
        j = 3
    if "545" in fichier :
        j = 4
    if "857" in fichier :
        j = 5

    print j 
    # map_gauss = fct.smooth(fichier, ma.sqrt(FWMH[0,1]**2 - FWMH[j,1]**2))
    #i = i + 1

