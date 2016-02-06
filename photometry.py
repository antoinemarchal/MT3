from photutils import aperture_photometry
from astropy.utils.data import download_file
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import function as fct

path ='patch_SZ/SZ/1_PSZ2 G000.13+78.04.fits'
PSZ = "PSZ2v1.fits"
NAME, GLON, GLAT = fct.coord_SZ(PSZ)
map_data = pyfits.getdata(path)
plt.figure()
plt.imshow(map_data)

