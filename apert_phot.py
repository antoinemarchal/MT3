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
from numpy import sqrt, log

def Gamma2sigma(Gamma):
    return Gamma / (sqrt(log(256)))

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

n_cluster = 1383 #Number of cluster selected
files = "patch_SZ/SZ/reduced_filenames.txt"
files = "patch_SZ/SZ/filenames.txt"
path  = "patch_SZ/SZ/"

threshold = 0.4

flux, redshift, MSZ, rcrit, rp = ap.do_photometry(
    n_cluster, files, path, threshold,0
)

slct_redshift = []
slct_flux     = []
slct_msz      = []
slct_rcrit    = []
slct_rp       = []

l = 0

pouet = 0
=======
sigma = Gamma2sigma(11.2466)
x = np.linspace(0,200,200)
gauss = gaussian(x, 0., sigma)
>>>>>>> b4a0ea0f4faa2fde9822e72554fbdcf0b7d78a2a
for i in range(len(redshift)):
    if redshift[i] >= 0. and MSZ[i] != 0. :
        slct_redshift.append(redshift[i])
        slct_flux.append(flux[i])
        slct_msz.append(MSZ[i])
        slct_rcrit.append(rcrit[i])
        slct_rp.append(rp[i])
        l +=1

        
            
        if rcrit[i] == 5 :
            fig_3   = plt.figure(figsize=(10,6))
            ax_3    = fig_3.add_subplot(1, 1, 1)
            ax_3.set_xlim(0, 35)
            print rcrit[i]
            plt.plot(rp[i])
            plt.plot(x, gauss)
            plt.show()
n_cl = l
print n_cl
k = 0

ar_rc = np.asarray(slct_rcrit)
rp = np.asarray(slct_rp)
rc = 0
nb_indexes = [0]

## determination du profil type en fonction du rayon critique

for k in range(np.max(slct_rcrit)) :
    indexes = np.where(ar_rc == k)
    indexes = np.asarray(indexes)
    n1,n2 = indexes.shape
    nb_indexes = np.append(nb_indexes,n2) 
    #retirer la premiere ligne a la sortie du bloc 
    if k == 0 :
        med_profile = np.zeros([1,150])
    
    if n2 == 0 :

        med_profile_temp = np.zeros([1,150]) 
        med_profile = np.vstack((med_profile,med_profile_temp))
        rc = rc + 1
        continue
    
    rp_rc = rp[indexes,:]
    rp_rc = rp_rc.reshape(n2,150)
    
    test_1d = 0
    for j in range(n2) :
        if j == 0 :
           matrix_rc = rp_rc[0,:]
           test_1d = 1  
        else:
            test_1d = 0
            matrix_rc = np.vstack((matrix_rc,rp_rc[j,:]))

    if test_1d == 1 :
        med_profile_temp = matrix_rc
    else : 
        med_profile_temp = np.median(matrix_rc, axis = 0)
        
 
    med_profile = np.vstack((med_profile,med_profile_temp))
    
# med_profile [i,:] designe le profile median  pour rc = i
# med_profile [0,:] = 0
# Si pas de rc a i , med_profile[i,:] = 0


#print nb_indexes
#a = 0.
#if a == 1 : 
#    for i in range(20):
#	if nb_indexes[i] >= 10  : 
#   	    plt.plot(med_profile[i,:])
#    plt.plot([0,150],[0.4,0.4])
#    plt.show()


"""---------------------------------------------------
---Save results
---------------------------------------------------"""
with open('results/flux.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(flux)
output.close()

with open('results/slct_redshift.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(slct_redshift)
output.close()
with open('results/slct_flux.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(slct_flux)
output.close()
with open('results/slct_msz.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(slct_msz)
output.close()
with open('results/slct_rcrit.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(slct_rcrit)
output.close()
with open('results/med_profile.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(med_profile)
output.close()
with open('results/nb_indexes.pkl', 'wb') as output:
    mon_pickler = pickle.Pickler(output)
    mon_pickler.dump(nb_indexes)
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

