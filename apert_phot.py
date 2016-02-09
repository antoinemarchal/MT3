import sys
sys.path.append('photutils')
from photutils import aperture_photometry
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import function as fct

path ='patch_SZ/SZ/23_PSZ2 G006.76+30.45.fits'
PSZ = "PSZ2v1.fits"
NAME, GLON, GLAT = fct.coord_SZ(PSZ)
map_data = pyfits.getdata(path)
plt.figure()
plt.imshow(map_data)

