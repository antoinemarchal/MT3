from astropy.io import fits as pyfits
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os

plt.ion()
filename = "maps_2015/HFI_SkyMap_143_2048_R2.00_full.fits"
filename_smooth = "maps_smooth/HFI_SkyMap_143_2048_R2.00_full.fits"
map,header_1 = hp.read_map(filename,h=True)
map_smooth,header_2 = hp.read_map(filename_smooth,h=True)
hp.mollview(map - map_smooth, norm='hist')
