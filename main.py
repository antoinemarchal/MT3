from astropy.io import fits as pyfits
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma

plt.ion()
map = hp.read_map('maps_2015/HFI_SkyMap_100_2048_R2.00_full.fits')
#map_size = hp.get_map_size(map)
#hp.mollview(map, norm='hist')
#print(hp.npix2nside(map_size))

fwhm = np.radians(9.66/60)
print fwhm
map_gauss = hp.smoothing(map, fwhm)
hp.mollview(map_gauss, norm='hist')
