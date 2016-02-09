from astropy.io import fits as pyfits
import astropy.table as pytabs
import matplotlib.pyplot as plt
import pyfits as pf
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

def plot_weight(GLAT, st_w, n_cluster):
    plt.figure()
    plt.suptitle('fixme', size=16)
    #bins = np.linspace(-0.1, 0.1, 100)
    plt.subplot(2,3,1)
    plt.hist(st_w[:n_cluster,0], facecolor='b')
    plt.title('100 GHz')

    plt.subplot(2,3,2)
    plt.hist(st_w[:n_cluster,1], facecolor='g')
    plt.title('143 GHz')

    plt.subplot(2,3,3)
    plt.hist(st_w[:n_cluster,2], facecolor='r')
    plt.title('217 GHz')

    plt.subplot(2,3,4)
    plt.hist(st_w[:n_cluster,3], facecolor='c')
    plt.title('353 GHz')

    plt.subplot(2,3,5)
    plt.hist(st_w[:n_cluster,4], facecolor='m')
    plt.title('545 GHz')
    
    plt.subplot(2,3,6)
    plt.hist(st_w[:n_cluster,5], facecolor='k')
    plt.title('857 GHz')
    
    """------------------------------"""
    plt.figure()
    plt.suptitle('weight as function of frequency fixme', size=16)
    plt.subplot(2,3,1)
    plt.plot(GLAT[:n_cluster],st_w[:n_cluster,0], "b.")
    plt.title('100 GHz')
    
    plt.subplot(2,3,2)
    plt.plot(GLAT[:n_cluster],st_w[:n_cluster,1], "g.")
    plt.title('143 GHz')
    
    plt.subplot(2,3,3)
    plt.plot(GLAT[:n_cluster],st_w[:n_cluster,2], "r.")
    plt.title('217 GHz')
    
    plt.subplot(2,3,4)
    plt.plot(GLAT[:n_cluster],st_w[:n_cluster,3], "c.")
    plt.title('353 GHz')
    
    plt.subplot(2,3,5)
    plt.plot(GLAT[:n_cluster],st_w[:n_cluster,4], "m.")
    plt.title('545 GHz')
    
    plt.subplot(2,3,6)
    plt.plot(GLAT[:n_cluster],st_w[:n_cluster,5], "k.")
    plt.title('857 GHz')
    return 0 
