from astropy.io import fits as pyfits
import astropy.table as pytabs
import matplotlib.pyplot as plt
import pyfits as pf
import numpy as np
import sys
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
    cat      = pyfits.getdata(filename)
    data     = pytabs.Table(cat)
    NAME     = data['NAME']
    GLON     = data['GLON']
    GLAT     = data['GLAT']
    REDSHIFT = data['REDSHIFT']
    MSZ      = data['MSZ']
    return (NAME,GLON,GLAT,REDSHIFT,MSZ)

def save_fits(name, redshift, masse, patch,indice):
    path = "patch_SZ/SZ/"
    filename = path + str(indice+1)+'_' + name + ".fits"
    if os.path.isfile(filename) == 1:
        os.remove(filename)
    a        = np.array([redshift])
    b        = np.array([masse])
    hdu      = pyfits.PrimaryHDU(patch.value)
    hdulist  = pyfits.HDUList([hdu])
    col_1    = pyfits.Column(name='REDSHIFT', format='E', array=a)
    col_2    = pyfits.Column(name='MSZ', format='E', array=b)
    cols     = pyfits.ColDefs([col_1,col_2])
    tbhdu    = pyfits.BinTableHDU.from_columns(cols)
    prihdr   = pyfits.Header()
    prihdr['COMMENT'] = "Commentary fixme"
    prihdu = pyfits.PrimaryHDU(header=prihdr)
    hdulist.append(tbhdu)
    hdulist.writeto(filename)
    return 0

def plot_map(NAME, TSZ_map, CMB_KSZ_map):
    #"""---------------Plot CMB+KSZ and TSZ--------------"""
    fig = plt.figure(figsize=(10,4))
    ax = fig.add_subplot(1, 1, 1)
    plt.suptitle(NAME, size=16)
    plt.subplot(1,2,1)
    imgplot_1 = plt.imshow(TSZ_map, origin='lower',
                           interpolation='none')
    plt.title('TSZ map')

    plt.subplot(1,2,2)
    imgplot_2 = plt.imshow(CMB_KSZ_map, origin='lower',
                           interpolation='none')
    plt.title('CMB+KSZ map')
    plt.savefig('results/sz_effect.eps')
    return fig

def plot_w_hist(GLAT_slct, st_w, n_slct, style):
    """--------------------------------------------"""
    with plt.style.context(style, after_reset=True):
        freq  = [103.416 ,144.903, 222.598, 355.218, 528.4, 776.582]
        color = ['b','g', 'r', 'c', 'm', 'k'] 
        fig   = plt.figure(figsize=(16,9))
        ax    = fig.add_subplot(1, 1, 1)
        plt.suptitle('fixme', size=16)
        for i in range(6):
            bins = np.linspace(np.min(st_w[:n_slct,i]),
                               np.max(st_w[:n_slct,i]), 20)
            plt.subplot(2,3,i+1)
            plt.hist(st_w[:n_slct,i], bins,histtype='step',color=color[i])
            plt.title(str(freq[i])+'GHz')
        plt.savefig('results/w_hist.pdf')
    plt.show()
    return 0

def plot_w_glat(GLAT_slct, GLAT_excl, st_w, st_w_excl,
                n_slct, n_excl, style):
    """--------------------------------------------"""
    #with plt.style.context(style, after_reset=True):
    freq  = [103.416 ,144.903, 222.598, 355.218, 528.4, 776.582]
    color = ['b.','g.', 'r.', 'c.', 'm.', 'k.'] 
    fig   = plt.figure(figsize=(16,9))
    ax    = fig.add_subplot(1, 1, 1)
    #plt.suptitle('weight as function of frequency fixme --'
    #             + str(n_slct) + 'Selected' + str(n_excl)
    #             + 'exclu', size=16)
    for i in range(6):
        plt.subplot(2,3,i+1)
        plt.plot(GLAT_slct[:n_slct], st_w[:n_slct,i], color[i],
                 GLAT_excl[:n_excl], st_w_excl[:n_excl,i], "y.")
        plt.title(str(freq[i])+'GHz')
        plt.xlabel('GLAT')
        plt.ylabel('Weight')
        plt.savefig('results/w_lat.pdf')
    plt.show()
    return 0

def  plot_w_full_patch(w_full, moy_w_slct, moy_w, std_w_slct
                       , std_w, W2, n_slct, style):
    #with plt.style.context(style, after_reset=True):
    freq = [103.416 ,144.903, 222.598, 355.218, 528.4, 776.582]    
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(1, 1, 1)
    plt.errorbar(freq, moy_w, marker="o", xerr=0., yerr=std_w,
                 color='cyan', label='1653 Cluster')
    plt.errorbar(freq, moy_w_slct, marker="o", xerr=0., yerr=std_w_slct,
                 color='blue', label=str(n_slct)+'Cluster')
    #plt.plot(freq, w_full, "r:o", lw=1, label='Full sky')
    plt.plot(freq, W2, "m:o", lw=1, label='Full sky') 
    plt.plot([0, 790], [0, 0], 'g--', lw=1)
    plt.legend()
    plt.xlabel('Frequence [GHz]')
    plt.ylabel('Weight')
    plt.grid()
    plt.savefig('results/w_full.pdf')
    plt.show()
    return 0

def progress(count, total, suffix=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s --- %s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()
