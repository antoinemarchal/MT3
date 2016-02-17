import astropy.table as ta
import astropy.io.fits as af 
cat, hdr = af.getdata("PSZ2v1.fits", header=True)
data = ta.Table(cat)
print data
print hdr
MSZ = data['MSZ']
print MSZ
