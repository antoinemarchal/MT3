from astropy.io import fits as pyfits
import astropy.table as pytabs
import matplotlib.pyplot as plt
import pyfits as pf
import numpy as np
import os


def coord_SZ(filename):
    """--------------------------------------------------
    --- coord_SZ : Getting back useful information of
                   galaxy cluster in SZ catalog 
    
                   Parameters : filename  = .fits file
    
                   Return     : NAME = ID of cluster  
                                GLON = galactic longitude
                                GLAT = galactic latitude
    -----------------------------------------------------"""
    cat  = pyfits.getdata(filename)
    data = pytabs.Table(cat)
    NAME = data['NAME']
    GLON = data['GLON']
    GLAT = data['GLAT']
    return (NAME,GLON,GLAT)

def save_fits(name, patch,indice):
    path = "patch_SZ/SZ/"
    filename = path + str(indice+1)+'_' + name + ".fits"
    if os.path.isfile(filename) == 1:
        os.remove(filename)    
    pf.writeto(filename, patch.value)
    return 0

def plot_map(NAME, TSZ_map, CMB_KSZ_map):
    #"""---------------Plot CMB+KSZ and TSZ--------------"""
    plt.suptitle(NAME, size=16)

    plt.subplot(1,2,1)
    plt.imshow(TSZ_map, origin='lower', interpolation='none')
    plt.title('TSZ map')
    
    plt.subplot(1,2,2)
    plt.imshow(CMB_KSZ_map, origin='lower', interpolation='none')
    plt.title('CMB+KSZ map')
    return 0

def plot_weight(GLAT_slct, GLAT_excl, st_w, st_w_excl, n_slct, n_excl):
    plt.figure()
    plt.suptitle('fixme', size=16)
    #bins = np.linspace(-0.1, 0.1, 100)
    plt.subplot(2,3,1)
    plt.hist(st_w[:n_slct,0], facecolor='b')
    plt.title('100 GHz')

    plt.subplot(2,3,2)
    plt.hist(st_w[:n_slct,1], facecolor='g')
    plt.title('143 GHz')

    plt.subplot(2,3,3)
    plt.hist(st_w[:n_slct,2], facecolor='r')
    plt.title('217 GHz')

    plt.subplot(2,3,4)
    plt.hist(st_w[:n_slct,3], facecolor='c')
    plt.title('353 GHz')

    plt.subplot(2,3,5)
    plt.hist(st_w[:n_slct,4], facecolor='m')
    plt.title('545 GHz')
    
    plt.subplot(2,3,6)
    plt.hist(st_w[:n_slct,5], facecolor='k')
    plt.title('857 GHz')
    
    """------------------------------"""
    plt.figure()
    plt.suptitle('weight as function of frequency fixme --'
                 + str(n_slct) + 'Selected' + str(n_excl)
                 + 'exclu', size=16)
    plt.subplot(2,3,1)
    plt.plot(GLAT_slct[:n_slct], st_w[:n_slct,0], "b.", GLAT_excl[:n_excl]
             ,st_w_excl[:n_excl,0], "y.")
    plt.title('100 GHz')
    
    plt.subplot(2,3,2)
    plt.plot(GLAT_slct[:n_slct], st_w[:n_slct,1], "g.", GLAT_excl[:n_excl]
             ,st_w_excl[:n_excl,1], "y.")
    plt.title('143 GHz')
    
    plt.subplot(2,3,3)
    plt.plot(GLAT_slct[:n_slct], st_w[:n_slct,2], "r.", GLAT_excl[:n_excl]
             ,st_w_excl[:n_excl,2], "y.")
    plt.title('217 GHz')
    
    plt.subplot(2,3,4)
    plt.plot(GLAT_slct[:n_slct], st_w[:n_slct,3], "c.", GLAT_excl[:n_excl]
             ,st_w_excl[:n_excl,3], "y.")
    plt.title('353 GHz')
    
    plt.subplot(2,3,5)
    plt.plot(GLAT_slct[:n_slct], st_w[:n_slct,4], "m.", GLAT_excl[:n_excl]
             ,st_w_excl[:n_excl,4], "y.")
    plt.title('545 GHz')
    
    plt.subplot(2,3,6)
    plt.plot(GLAT_slct[:n_slct],st_w[:n_slct,5], "k.", GLAT_excl[:n_excl]
             ,st_w_excl[:n_excl,5], "y.")
    plt.title('857 GHz')
    return 0 

def  plot_w_full_patch(w_full, moy_w_slct, moy_w, std_w_slct, std_w, W2):
    freq = [103.416 ,144.903, 222.598, 355.218, 528.4, 776.582]    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.errorbar(freq, moy_w, marker="o", xerr=0., yerr=std_w,
                 color='cyan', label='1653 Cluster')
    plt.errorbar(freq, moy_w_slct, marker="o", xerr=0., yerr=std_w_slct,
                 color='blue', label='1360 Cluster')
    plt.plot(freq, w_full, "r:o", lw=1, label='Full sky')
    plt.plot(freq, W2, "m:o", lw=1, label='Full sky perso') 
    plt.plot([0, 790], [0, 0], 'g--', lw=1)
    plt.legend()
    plt.show
    return 0
