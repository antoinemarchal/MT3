import sys
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import in_output as inout
import numpy.ma as ma
import math as ma 



def sector_mask(shape,centre,radius,angle_range):
    
    #retourne un masque booleen sur un secteur circulaire
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
        #plt.show()
        plt.subplot(2,2,2)
        plt.imshow(data_circle,origin = 'lower',interpolation = 'none')
        #plt.show()
        plt.subplot(2,2,3)
        plt.imshow(data_ring,origin = 'lower',interpolation = 'none')
        #plt.show()
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
    
    #convert r to int ==> definit la taille du bin
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    
    #Normalisation
    radialprofile = (radialprofile - np.min(radialprofile)) \
                    / (np.max(radialprofile) - np.min(radialprofile))
                           
    #### Calcul d'un rayon caracteristique rc :
    #### 90 % du flux en partant du centre
    r_90 = np.where(radialprofile >= 0.1) #FIXME
    r_90 = np.asarray(r_90)
    r_90 = np.ravel(r_90)
    for i in range(len(r_90)-1):
        if r_90[i+1] != r_90[i]+1:
            rc = r_90[i]
            break
        else:
            rc = np.max(r_90)
            
    if plot == 1 : 
        fig = plt.figure(4)
        ax = fig.add_subplot(1, 1, 1)
        plt.plot(radialprofile)
        plt.plot([rc, rc], [0,1], 'r--', lw=2)
        plt.plot([0, 70], [0.1, 0.1], 'g--', lw=2)
        ax.axvspan(0, rc, alpha=0.5, color='grey')
        plt.xlabel('Radius [pix]')
        plt.ylabel('Flux')
        #plt.show()
    return (radialprofile,rc) 

#######################################################

def get_flux(data_circle,data_ring):
#retourne la valeur du flux de la cource en tenant compte 
#du bruit moyen dans l'anneau
	raw_flux =  np.sum (data_circle)

	pixel_ring =np. where(data_ring > 0 )
	pixel_ring = np.ravel(pixel_ring)
	npixel_ring = len(pixel_ring)
	
	pixel_circle = np.where(data_circle > 0 )
	pixel_circle = np.ravel(pixel_circle)
	npixel_circle = len(pixel_circle)

	avg_bckd = np.sum(data_ring)/npixel_ring
	
	flux = raw_flux - npixel_circle*avg_bckd

	return flux,npixel_circle,npixel_ring
#######################################################

def do_photometry():
	
	filenames =  open ("patch_SZ/SZ/filenames.txt")
	path = "patch_SZ/SZ/"
        k = 0
	for line in filenames:
            inout.progress(k, 1320, 'Cluster')
	    patch = path + line.strip()
	    data = pyfits.getdata(patch)
	    n1,n2 = data.shape
	    centre = (n1/2,n2/2)
	    profile,rc = radial_profile(data,centre,0)
	    #print rc
	    ##FIXME
	    ##modifier les valeur de rayon pour les anneaux
	    ### !!!! ne jamais metre 1 en dernier argument####
	    data_circle,data_ring = phot_mask(data,rc,35,45,0)
	    flux =  get_flux(data_circle,data_ring)
            k +=1
	return 0 




patch ='patch_SZ/SZ/1441_PSZ2 G305.76+44.79.fits'
#patch = 'patch_SZ/SZ/121_PSZ2 G033.97-76.61.fits'
#patch = 'patch_SZ/SZ/276_PSZ2 G066.26+20.82.fits'
#patch = 'patch_SZ/SZ/348_PSZ2 G081.22-41.95.fits'
#patch = 'patch_SZ/SZ/1592_PSZ2 G340.94+35.07.fits'
data = pyfits.getdata(patch)
n1,n2 = data.shape
centre = (n1/2,n2/2)

profile,rc = radial_profile(data,centre,1)
print "rayon crit",rc
data_circle,data_ring = phot_mask(data,rc ,35,45,1)

flux = get_flux(data_circle,data_ring)

a = do_photometry()


#frac_in,frac_out = lobe_frac(data,data_circle)
#print "fracction du lobe dedans et dehors",frac_in, frac_out
#area_circle,area_ring = area(data_circle,data_ring)
#print "surfaces cercle et anneau : ", area_circle,area_ring



##################################################################################################
# penser a mettre des valeurs aberentes qunad : divisionpar 0 ou mauvaise valeur double scalaire 
