from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

import fonction as fct

plt.ion()
filename_smooth = "maps_smooth/HFI_SkyMap_143_2048_R2.00_full.fits"
map_smooth,header = hp.read_map(filename_smooth,h=True)
#hp.mollview(map_smooth, norm='hist')

"""
Creation of wcs objects for patch

FIXME information

"""
PSZ = "PSZ2v1.fits"
RA,DEC = fct.coord_sz(PSZ)

w = wcs.WCS(naxis=2)
naxis1 = 2
naxis2 = 3
w.wcs.crpix = [RA[0],DEC[0]]
# taille d'un pixel dans le ciel
w.wcs.cdelt = numpy.array([-0.066667, 0.066667])
#
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

#print RA
#(x,y) = w.wcs_world2pix(RA,DEC,0) # back to pixel coordinate
#print x
