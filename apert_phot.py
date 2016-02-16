import sys
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import in_output as inout
import numpy.ma as ma
import math as ma
import astropy.table as pytabs



def sector_mask(shape,centre,radius,angle_range):
    
    #retourne un masque booleen sur un secteur circulaire
    #shape      = data.shape
    #centre      = tuple
    #angle_range = tuple
    
    x,y = np.ogrid[:shape[0],:shape[1]]
    cx,cy = centre
    tmin,tmax = np.deg2rad(angle_range)

    # mensure stop angle > start angle
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
        plt.figure(1)
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
	#la normalisation me parait bidon
    radialprofile = (radialprofile - np.min(radialprofile)) \
                    / (np.max(radialprofile) - np.min(radialprofile))
                           
    #### Calcul d'un rayon caracteristique rc :
    #### 60 % du flux en partant du centre
	##valeur seuil a discuter...
    threshold = 0.4
    r_60 = np.where(radialprofile >= threshold)
    r_60 = np.asarray(r_60)
    r_60 = np.ravel(r_60)
    for i in range(len(r_60)-1):
        if r_60[i+1] != r_60[i]+1:
            rc = r_60[i]
            break
        else:
	#discutable si profile non monotone
            rc = np.max(r_60)

    if plot == 1 :
        style = 'grayscale'
        with plt.style.context(style, after_reset=True):
            fig = plt.figure(2)
            ax = fig.add_subplot(1, 1, 1)
            plt.plot(radialprofile, 'b')
            plt.plot([rc, rc], [0,1], 'r--', lw=2)
            plt.plot([0, 70], [threshold, threshold], 'g--', lw=2)
            ax.axvspan(0, rc, alpha=0.4, color='grey')
            plt.xlabel('Radius [pix]')
            plt.ylabel('Flux')
    #plt.show()
    return (radialprofile,rc) 

#######################################################

def get_flux(data_circle,data_ring):
#retourne la valeur du flux de la cource en tenant compte 
#du bruit moyen dans l'anneau
    raw_flux =  np.sum (data_circle)

    pixel_ring =np.where(data_ring != 0 )
    pixel_ring = np.ravel(pixel_ring)
    npixel_ring = len(pixel_ring)
   
    
    pixel_circle = np.where(data_circle != 0 )
    pixel_circle = np.ravel(pixel_circle)
    npixel_circle = len(pixel_circle)

    avg_bckd = np.sum(data_ring)/npixel_ring
    
    flux = raw_flux - npixel_circle*avg_bckd
    
    return flux#,npixel_circle,npixel_ring
#######################################################

def do_photometry(n_cluster):
    filenames =  open ("patch_SZ/SZ/filenames.txt")
    path = "patch_SZ/SZ/"
    k = 0
    flux      = []
    redshift  = []
    for line in filenames:
        inout.progress(k, n_cluster, 'Cluster')
	patch      = path + line.strip()
	data       = pyfits.getdata(patch)
        cat        = pyfits.getdata(patch,1)
        hdr        = pytabs.Table(cat)
        rd         = hdr['REDSHIFT']
    	n1,n2      = data.shape
	centre     = (n1/2,n2/2)
	profile,rc = radial_profile(data,centre,0)
       
	##FIXME
	##modifier les valeur de rayon pour les anneaux
	### !!!! ne jamais metre 1 en dernier argument####
	data_circle,data_ring = phot_mask(data,rc,35,45,0)
        pouet = get_flux(data_circle,data_ring)
        #if pouet >= 0.03 :
        #    plt.figure()
        #    plt.imshow(data_circle)
        #    plt.show()
        
        #if rc <=13 : 
        flux.append(get_flux(data_circle,data_ring))
        redshift.append(rd[0])
        k += 1
    return flux, redshift
    

#**************************************************************************
#**************************************************************************
#**************************************************************************
#**************************************************************************


PSZ = "PSZ2v1.fits"
NAME,GLON,GLAT, REDSHIFT = inout.coord_SZ(PSZ)

n_cluster = 1391  
flux, redshift = do_photometry(n_cluster)
l = 0
rslt = np.zeros((len(flux),len(flux)))
for i in range(len(redshift)):
    if redshift[i] >= 0.:
        rslt[l][0]= redshift[i]
        rslt[l][1]= flux[i]
        l +=1
n_cl = l

moy   = np.mean(rslt[:,1])
std = np.std(rslt[:,1])
out_rslt = np.zeros((n_cl,n_cl))
in_rslt  = np.zeros((n_cl,n_cl))
j = 0
l = 0
for i in range(n_cl):
    RD = np.absolute((rslt[i][1]-moy) / std)
    if RD > 1.:
        out_rslt[j][0]=rslt[i][0]
        out_rslt[j][1]=rslt[i][1]
        j += 1
    else:
        in_rslt[l][0]=rslt[i][0]
        in_rslt[l][1]=rslt[i][1]
        l += 1
n_in  = l
n_out = j
print '-'
print str(n_in)  + 'Cluster selected'
print str(n_out) + 'Cluster excluded'

"""----------Plot/Results Study flow/redshift------------"""
bins_r = np.linspace(0., np.max(in_rslt[:,0]), 40)
bins_f = np.linspace(np.min(in_rslt[:,1]),
                     np.max(in_rslt[:,1]), 40) #FIXME problem Nan
color = ['b.','g.', 'r.', 'c.', 'm.', 'k.', 'y.'] 
fig   = plt.figure()
ax    = fig.add_subplot(1, 1, 1)
plt.subplot(2,2,1)
plt.xlabel('z')
plt.ylabel('Flux')
plt.plot(in_rslt[:n_in,0], in_rslt[:n_in,1], color[2],
         out_rslt[:n_out,0], out_rslt[:n_out,1], color[6])
plt.subplot(2,2,2)
plt.xlabel('z')
plt.ylabel('N cluster')
plt.hist(in_rslt[:n_in,0], bins_r, facecolor='b')
plt.subplot(2,2,3)
plt.xlabel('Flux')
plt.ylabel('N cluster')
plt.hist(in_rslt[:n_in,1], bins_f, facecolor='g')
plt.plot([moy, moy], [0,85], 'r--', lw=2)
plt.savefig('results/study_flow.png')
plt.show()



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

