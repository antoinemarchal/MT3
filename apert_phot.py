import sys
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import function as fct
import numpy.ma as ma
import math as ma 



def sector_mask(shape,centre,radius,angle_range):
    
    #retourne un masque booleen sur un secteur sirculaire
    #shape      = data.shape
    #centre      = tuple
    #angle_range = tuple
    
    x,y = np.ogrid[:shape[0],:shape[1]]
    cx,cy = centre
    tmin,tmax = np.deg2rad(angle_range)

    # ensure stop angle > start angle
    if tmax < tmin:
            tmax += 2*np.pi

    # convert cartesian --> polar coordinates
    r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy)
    theta = np.arctan2(x-cx,y-cy) - tmin

    # wrap angles between 0 and 2*pi
    theta %= (2*np.pi)

    # circular mask
    circmask = r2 <= radius*radius

    # angular mask
    anglemask = theta <= (tmax-tmin)

    return circmask*anglemask
#####################################################################

def phot_mask(data,r_circle,rin,rout,plot):
# data     = donnees extraites du patch
#r_circle  = rayon du cercle
#rin       = rayon interne de l'anneau
#rout      = rayon externe de l'anneau
#plot      = 1 ou 0 si =o on trace rien

   
    data_f = data
    data_circle = np.copy(data)
    data_rin = np.copy(data)
    data_rout = np.copy(data)
    
    centre_x,centre_y = data.shape
    x = centre_x/2
    y = centre_y/2
    
    mask_circle = sector_mask(data_circle.shape,(x,y),r_circle,(0,360))
    mask_rin    = sector_mask(data_rin.shape,(x,y),rin,(0,360))
    mask_rout   = sector_mask(data_rout.shape,(x,y),rout,(0,360))

    data_circle[~mask_circle]= 0
    data_rin[~mask_rin]      = 0 
    data_rout[~mask_rout]    = 0

    data_ring =  data_rout - data_rin
    
    if plot == 1: 
        plt.figure(2)
        plt.subplot(2,2,1)
        plt.imshow(data,origin = 'lower',interpolation = 'none')
        plt.show()
        plt.subplot(2,2,2)
        plt.imshow(data_circle,origin = 'lower',interpolation = 'none')
        plt.show()
        plt.subplot(2,2,3)
        plt.imshow(data_ring,origin = 'lower',interpolation = 'none')
        plt.show()
        plt.subplot(2,2,4)
        plt.imshow(data_circle + data_ring,origin='lower',interpolation='none')
        plt.show()
    
    
    return (data_circle,data_ring)
##################################################################

def lobe_frac(data,data_circle):
    # retourne la fraction du signal SZ contenu dans et a l'exterieure
    # de rayon r_in
    #data          = patch complet 
    #data_circle   = ouverture photometrique
    #data_ring     = anneau autour de l'ouverture photometrique

    tot =np.sum(data)
    int = np.sum(data_circle)

    frac_in = int/tot
    frac_out = 1-frac_in
    
    return (frac_in,frac_out) 
##################################################################

def area(data_circle, data_ring):
    #retourne les valeurs des surfaces en pixels a l'interieure de
    # l'anneau et de l'ouverture circulaire.

    area_circle = np.count_nonzero(data_circle)
    area_ring   = np.count_nonzero(data_ring)
    
    return (area_circle,area_ring)

###################################################################


def radial_profile(data, center,plot):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    
    #convert r to int ==> definit la tille du bin
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr

    if plot ==1 : 
        plt.figure(4)
        plt.plot(radialprofile)
        plt.show()

    #### calcul d'un rayon caracteristique rc :
    max = np.max(radialprofile)
    mean_noise = np.mean(radialprofile[40:50]) #FIXME
  
    r2 = np.where(radialprofile-mean_noise >= 0.3*max)
    r2 = np.asarray(r2)
    r2 = np.ravel(r2)
    rc = r2[-1]

    return (radialprofile,rc) 



patch ='patch_SZ/SZ/23_PSZ2 G006.76+30.45.fits'
data = pyfits.getdata(patch)
n1,n2 = data.shape
centre = (n1/2,n2/2)

profile,rc = radial_profile(data,centre,1)
print "rayon crit",rc
data_circle,data_ring = phot_mask(data,rc ,35,45,1)
#### jusque la c'est bon


frac_in,frac_out = lobe_frac(data,data_circle)
print "fracction du lobe dedans et dehors",frac_in, frac_out
area_circle,area_ring = area(data_circle,data_ring)
print "surfaces cercle et anneau : ", area_circle,area_ring




