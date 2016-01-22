from astropy.io import fits as pyfits
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

plt.ion()
map = hp.read_map('maps_2015/HFI_SkyMap_100_2048_R2.00_full.fits')
hp.mollview(map, norm='hist')

