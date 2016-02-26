import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.font_manager import FontProperties
import numpy as np
from numpy import sqrt, log
import math as ma
import pickle
import plot 

def Gamma2sigma(Gamma):
    return Gamma / ( sqrt(log(256)))

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

flux            = pickle.load(open("results/flux.pkl","rb"))
slct_redshift   = pickle.load(open("results/slct_redshift.pkl","rb"))
slct_flux       = pickle.load(open("results/slct_flux.pkl","rb"))
slct_msz        = pickle.load(open("results/slct_msz.pkl","rb"))
slct_rcrit      = pickle.load(open("results/slct_rcrit.pkl","rb"))
med_profile     = pickle.load(open("results/med_profile.pkl","rb"))
nb_indexes      = pickle.load(open("results/nb_indexes.pkl","rb"))

moy   = np.mean(slct_flux)
std   = np.std(slct_flux)
moy_msz=np.mean(slct_msz)
"""----------Plot/Results Study flow/redshift------------"""
bins_r = np.linspace(0., 1., 100)
bins_f = np.linspace(0., 0.025, 60)
bins_m = np.linspace(np.min(slct_msz),
                     np.max(slct_msz), 40)
bins_rc = np.linspace(np.min(slct_rcrit),
                      np.max(slct_rcrit), 32)
fontP = FontProperties()
fontP.set_size('xx-small')
color = ['b.','g.', 'r.', 'c.', 'm.', 'k.', 'y.'] 
fig_1   = plt.figure(figsize=(12,9))
ax_1    = fig_1.add_subplot(1, 1, 1)
plt.locator_params(nbins=4)
plt.subplot(2,2,1)
plt.xlabel('$z$')
plt.ylabel('$Flux$')
plt.plot(slct_redshift, slct_flux, color[5])

plt.subplot(2,2,2)
plt.xlabel('$z$')
plt.ylabel('$R_c$')
plt.plot(slct_redshift, slct_rcrit, color[5])

plt.subplot(2,2,3)
plt.xlabel('$z$')
plt.ylabel('$N_{cluster}$')
plt.hist(slct_redshift, bins_r,
            histtype='step', color='b')
plt.tight_layout()
plt.subplot(2,2,4)
plt.xlabel('$Flux$')
plt.ylabel('$N_{cluster}$')
plt.hist(flux, bins_f,
         histtype='step', color='b', label='$Unknow$ $z$')
plt.hist(slct_flux, bins_f,
         histtype='step', color='g', label='$Know$ $z$')
plt.plot([moy, moy], [0,130], 'r--', lw=2,
         label=r'$\bar F$ $for$ $know$ $z$')
plt.legend(loc=1, numpoints=1)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='xx-small')
plt.savefig('results/rslt_1.pdf', format='pdf')

fig_2   = plt.figure(figsize=(12,9))
ax_2    = fig_2.add_subplot(1, 1, 1)
plt.subplot(2,2,1)
plt.xlabel('$M/10^{14} M \odot$')
plt.ylabel('$Flux$')
plt.plot(slct_msz, slct_flux, color[5])

plt.subplot(2,2,2)
plt.xlabel('$M/10^{14} M \odot$')
plt.ylabel('$N_{cluster}$')
plt.hist(slct_msz, bins_m,
            histtype='step', color='c')
plt.plot([moy_msz,moy_msz], [0,100], 'r--', lw=2)
plt.tight_layout()

plt.subplot(2,2,4)
plt.xlabel('$R_c$')
plt.ylabel('$N$')
plt.hist(slct_rcrit, bins_rc,
            histtype='step', color='g')
plt.tight_layout()
plt.savefig('results/rslt_2.pdf', format='pdf')

fig_3   = plt.figure(figsize=(10,6))
ax_3    = fig_3.add_subplot(1, 1, 1)
ax_3.set_xlim(0, 35)
height = 1
mean = 0
sigma = Gamma2sigma(11.2466)
x = np.linspace(0,40,100)
gauss = gaussian(x, 0., sigma)
plt.plot(x,gauss, 'k--', label='PSF')
for i in range(20):
    if nb_indexes[i] >= 2 :
   	plt.plot(med_profile[i,:], label="$R_{c} = $"+str(i))
plt.plot([0,150],[0.4,0.4])
plt.plot([0,150],[0.5,0.5])
plt.xlabel('$R$ $[pix]$')
plt.ylabel('$Normalized$ $Flux$')
plt.legend(loc=1, numpoints=1)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='small')
plt.savefig('results/radius.pdf', format='pdf')

plt.show()
