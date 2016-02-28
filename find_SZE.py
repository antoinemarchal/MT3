import argparse
from astropy.io import fits as pyfits
from astropy import wcs 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import pickle
import sys
import os

import mod_ilc as ilc
import in_output as inout

plt.ion()

parser = argparse.ArgumentParser()
parser.add_argument("n_cluster", help="Choose the number of cluster",
                    type=int)
parser.add_argument("-p", "--plot", help="Plot results",
                    action="store_true")
args = parser.parse_args()
if args.n_cluster > 1653:
    print 'There are only 1653 Clusters'
    sys.exit(1)

n_obs          = 6               #Number of HFI maps Planck
n_cluster      = args.n_cluster  #Number of cluster in catalog
patch_size     = 250

PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT, REDSHIFT, MSZ = inout.coord_SZ(PSZ)

st_w          = np.zeros((n_cluster,n_obs)) #weight of selected patch
st_w_excl     = np.zeros((n_cluster,n_obs)) #weight of excluded patch 
GLAT_slct     = np.zeros(n_cluster)         
GLAT_excl     = np.zeros(n_cluster)
fwhm          = np.loadtxt("fwhm_HFI.txt")
freq          = []
for i in range(n_obs):
    freq.append(fwhm[i,0])

"""-----------------------------------------------------
   Read the average value and standard deviation of 
   weight (compute and save one time) for each frequency 
   using weight.dat file  
--------------------------------------------------------"""
w = pickle.load(open("weight.dat","rb"))
moy_w      = []
std_w      = []
for i in range(n_obs):
    moy_w.append(np.mean(w[:n_cluster,i]))
    std_w.append(np.std(w[:n_cluster,i]))    
"""-----------------------------------------------------
   Reading the 6 smooth HFI maps Planck 
--------------------------------------------------------"""
unit_1 = open("files_HFI_full.txt")
origin = unit_1.tell()
path_1 = "maps_smooth/"
map_smooth = []
header     = []
for line in unit_1:
    sys.stdout = open(os.devnull, "w")
    filename_smooth = line.strip()
    map_smooth.append(
        hp.read_map(path_1 + filename_smooth)
    )
    sys.stdout = sys.__stdout__
unit_1.seek(origin)
"""-----------------------------------------------------
   Separation of SZ and CMB componant of 
   full sky using ilc method
   -----------------------------------------------------"""
mask = hp.read_map('galacticMask.fits')
map_mask = []
for i in range(6) :
    map_mask.append(map_smooth[i] * mask)
    
print 'ILC on full sky'
f_nu          = ilc.dist_SZ(freq)
a             = np.ones(n_obs)
a_t           = np.transpose(a)
b             = np.ones(n_obs)
b             = f_nu
b_t           = np.transpose(b)
CMB   = np.zeros(hp.get_map_size(map_smooth[0]))
TSZ   = np.zeros(hp.get_map_size(map_smooth[0]))
J = np.cov((map_smooth[0], map_smooth[1], map_smooth[2]
            , map_smooth[3], map_smooth[4], map_smooth[5]))
K = np.linalg.inv(J)
WK, WT = ilc.weight(K, a, a_t, b ,b_t)

# test with mask
L = np.cov((map_mask[0], map_mask[1], map_mask[2]
            , map_mask[3], map_mask[4], map_mask[5]))
P = np.linalg.inv(L)
WK_mask, WT_mask = ilc.weight(P, a, a_t, b ,b_t)

for i in range(n_obs) :
    CMB = CMB + (WK[i] * map_smooth[i])
    TSZ = TSZ + (WT[i] * map_smooth[i])
#FIXME we have not the same results with Simon
#and Alex approch
"""----------------------------------------------------"""
"""--------------------------------------------
   Creating patches on each cluster of catalog
   and extraction of SZ and CMB componant using 
   ILC method.
   --------------------------------------------"""
print "Creating patches of SZ effect : "
filename = "patch_SZ/SZ/filenames.txt"
if os.path.isfile(filename) == 1:
    os.remove(filename)
fichier = open(filename, "a")

moy_w_slct = []
std_w_slct = []
j = 0
l = 0
for k in range(n_cluster):
    inout.progress(k, n_cluster, 'Cluster selected = '+str(l+1))
    #Compute results of separation and weight for
    #each cluster.
    sys.stdout = open(os.devnull, "w")
    CMB_KSZ_map, TSZ_map, w_t, w_k = ilc.separation(
        k, map_smooth, patch_size, NAME, GLON,
        GLAT, freq
    )
    sys.stdout = sys.__stdout__
    
    """------Separation of cluster population------
       --------------------------------------------"""
    #Compute of relative deviation using standard
    #deviation (SD) and average value of weight.
    RD = []
    for i in range(n_obs-2):
        RD.append(
            np.absolute((w_t[i]-moy_w[i]) / std_w[i])
        )
    if args.plot:
        if k == 22:
            fig_map = inout.plot_map(NAME[k], TSZ_map, CMB_KSZ_map)

    #Condition with SD on each freq to select cluster
    #of our catalog.
    if (
            (RD[0] > 1.)| (RD[1] > 1.)|
            (RD[2] > 1.)| (RD[3] > 1.)
    ):
        GLAT_excl[j]   = GLAT[k]
        st_w_excl[j,:] = w_t
        j += 1
    else:
        GLAT_slct[l]   = GLAT[k]
        st_w[l,:]      = w_t
        inout.save_fits(NAME[k], REDSHIFT[k], MSZ[k], TSZ_map, k)
        fichier.write(str(k+1)+'_'
                      + NAME[k] + ".fits"+"\n")
        l += 1
n_slct = l
n_excl = j
fichier.close()
#with open('weight.dat', 'wb') as output:
#    mon_pickler = pickle.Pickler(output)
#    mon_pickler.dump(st_w)
"""---------------------------------------------------------"""
#Compute the average and SD of weight for selected cluster
#for each frequences of HFI maps Planck
for i in range(n_obs):
    moy_w_slct.append(np.mean(st_w[:n_slct,i]))
    std_w_slct.append(np.std(st_w[:n_slct,i]))

w_full_sky = np.loadtxt("w_full_sky.d")

if args.plot:
    print "Results :"
    #Plot of weight histogram and GLAT correlation
    #Plot of comparison with average weight of other group
    inout.plot_w_hist(
        GLAT_slct, st_w, n_slct, 'bmh'
    )
    inout.plot_w_glat(
        GLAT_slct, GLAT_excl, st_w, st_w_excl,
        n_slct, n_excl, 'ggplot'
        )
    inout.plot_w_full_patch(
        w_full_sky[:,0], moy_w_slct, moy_w, std_w_slct,
        std_w, WT, WT_mask, n_slct, 'ggplot'
    )
    




    
