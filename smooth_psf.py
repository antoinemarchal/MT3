from astropy.io import fits as pyfits
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os
import smooth_psf as fct

def smooth(filename, fwhm):
    """-------------------------------------------------
    --- smooth : Convoluted form by a gaussian.
    
                 Parameters : filename  = adress of map
    
                 Return     : map_gauss = convoluted map
                                          header
    ----------------------------------------------------"""
    map,header = hp.read_map(filename,h=True)
    fwhm_rad = np.radians(fwhm/60.)
    map_gauss = hp.smoothing(map, fwhm_rad)
    return (map_gauss, header)


"""------------------------------------------------------------
   ---Program : Convolution using smoothing function of healpy 
                to degraded HFI images of Planck before apply
                linear combinaison to study SZ effect
   ------------------------------------------------------------"""
   
unit_1 = open("filenames_HFI.txt")
FWMH = np.loadtxt("fwhm_HFI.txt")
i = 0
for line in unit_1 :
    fichier = line.strip()
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
        
    map_smooth, header = fct.smooth(
        fichier, ma.sqrt(FWMH[0,1]**2 - FWMH[j,1]**2)
    )
    
    if "545" in fichier :
        map_smooth /= 58.04 # MJy.s-1 to Tcmb convertion FIXME unit 
    if "857" in fichier : 
        map_smooth /= 2.27  # MJy.s-1 to Tcmb convertion FIXME unit

    i = i + 1
    path_1 = "maps_smooth/"
    if os.path.isfile(path_1 + fichier[10:]) == 1 :
        os.remove(path_1 + fichier[10:])
    hp.write_map(
        "maps_smooth/" + fichier[10:],map_smooth, extra_header=(header)
    )
        
for line in unit_1 :
    fichier = line.strip()
    if "100" in fichier :
        map,header = hp.read_map(fichier,h=True)
        hp.write_map(
            "maps_smooth/" + fichier[10:],map, extra_header=(header)
        )

