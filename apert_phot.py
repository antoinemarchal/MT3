import sys
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import in_output as inout
import numpy.ma as ma
import math as ma
import astropy.table as pytabs
import mod_ap as ap
import pickle

n_cluster = 1391
files = "patch_SZ/SZ/filenames.txt"
path  = "patch_SZ/SZ/"

r_in      = 35 #FIXME between rc and shape image
r_out     = 45
threshold = 0.4

flux, redshift, MSZ = ap.do_photometry(
    n_cluster, files, path, r_in, r_out, threshold
)

sort_redshift = []
sort_flux     = []
sort_msz      = []

l = 0
for i in range(len(redshift)):
    if redshift[i] >= 0. and MSZ[i] != 0. :
        sort_redshift.append(redshift[i])
        sort_flux.append(flux[i])
        sort_msz.append(MSZ[i])
        l +=1
n_cl = l
#FIXME look mass/redshift unknow ident?

moy   = np.mean(sort_flux)
std   = np.std(sort_flux)

RD_in_redshift  = []
RD_out_redshift = []
RD_in_flux      = []
RD_out_flux     = []
RD_in_msz       = []
RD_out_msz      = []

j = 0
l = 0
for i in range(n_cl):
    RD = np.absolute((sort_flux[i]-moy) / std)
    if RD > 2.:
        RD_out_redshift.append(sort_redshift[i])
        RD_out_flux.append(sort_flux[i])
        RD_out_msz.append(sort_msz[i])
        j += 1
    else:
        RD_in_redshift.append(sort_redshift[i])
        RD_in_flux.append(sort_flux[i])
        RD_in_msz.append(sort_msz[i])
        l += 1

n_in  = l
n_out = j
print '-'
print str(n_in)  + 'Cluster selected'
print str(n_out) + 'Cluster excluded'

"""---------------------------------------------------
---Save results
---------------------------------------------------"""
with open('results/sort_redshift.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(sort_redshift)
output.close()
with open('results/sort_flux.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(sort_flux)
output.close()
with open('results/sort_msz.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(sort_msz)
output.close()

with open('results/RD_in_redshift.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(RD_in_redshift)
output.close()
with open('results/RD_in_flux.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(RD_in_flux)
output.close()
with open('results/RD_in_msz.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(RD_in_msz)
output.close()

with open('results/RD_out_redshift.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(RD_out_redshift)
output.close()
with open('results/RD_out_flux.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(RD_out_flux)
output.close()
with open('results/RD_out_msz.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(RD_out_msz)
output.close()

##################################################################################################
# sources SZ non centrees ==> cercle photometrique bcp trop grand
# car le profile radial est fausse

# avec des techniques de roumains on peut voir que c'est probablement
# la bonne explication ==> bcp moin de source a tres grand flux 

#soit on degage les sources non centrees soit on fait un profil
# radial a partir du max (qu'on suppose etre le centre de la source)
# de la carte


## effet de resolution de source ? seul les amas avec un redshift faible on
# parfois un flux tres important.



#separer les populations et comparer a un catalogue de super amas 

# travailler sur des patch en 200*200  ==> penser a regarder l'influence de 
# tout ce bordel sur les poids

