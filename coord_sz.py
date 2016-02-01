import astropy.table as ta
import astropy.io.fits as af 
cat = af.getdata("PSZ2v1.fits")
data = ta.Table(cat)
print data
