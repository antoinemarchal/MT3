from astropy.io import fits as pf
import numpy as np

filename = "maps_2015/HFI_SkyMap_100_2048_R2.00_full.fits"
image = pf.open(filename)
var = image[0].data
print var
#pf.info(filename) 
#image, header = pf.getdata(filename,'image',header=True)

