import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from matplotlib.font_manager import FontProperties
import numpy as np
import pickle

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
bins_r = np.linspace(0., np.max(slct_redshift), 40)
bins_f = np.linspace(np.min(slct_flux),
                     np.max(slct_flux), 40)
bins_m = np.linspace(np.min(slct_msz),
                     np.max(slct_msz), 40)
bins_rc = np.linspace(np.min(slct_rcrit),
                      np.max(slct_rcrit), 32)
color = ['b.','g.', 'r.', 'c.', 'm.', 'k.', 'y.'] 
fig_1   = plt.figure(figsize=(10,6))
ax_1    = fig_1.add_subplot(1, 1, 1)
plt.locator_params(nbins=4)
plt.subplot(2,2,1)
plt.xlabel('$z$')
plt.ylabel('$Flux$')
plt.plot(slct_redshift, slct_flux, color[2])

plt.subplot(2,2,3)
plt.xlabel('$Flux$')
plt.ylabel('$N_{cluster}$')
plt.hist(slct_flux, bins_f,
            histtype='step', color='g')
plt.plot([moy, moy], [0,900], 'r--', lw=2)

plt.subplot(2,2,4)
plt.xlabel('$z$')
plt.ylabel('$N_{cluster}$')
plt.hist(slct_redshift, bins_r,
            histtype='step', color='b')
plt.tight_layout()
plt.savefig('results/flux_z.pdf', format='pdf')

fig_2   = plt.figure(figsize=(10,6))
ax_2    = fig_2.add_subplot(1, 1, 1)
plt.subplot(2,2,1)
plt.xlabel('$M/10^{14} M \odot$')
plt.ylabel('$Flux$')
plt.plot(slct_msz, slct_flux, color[4])

plt.subplot(2,2,2)
plt.xlabel('$M/10^{14} M \odot$')
plt.ylabel('$N_{cluster}$')
plt.hist(slct_msz, bins_m,
            histtype='step', color='c')
plt.plot([moy_msz,moy_msz], [0,100], 'r--', lw=2)
plt.tight_layout()
plt.savefig('results/flux_M.pdf', format='pdf')

plt.subplot(2,2,3)
plt.xlabel('$z$')
plt.ylabel('$R_c$')
plt.plot(slct_redshift, slct_rcrit, color[0])

plt.subplot(2,2,4)
plt.xlabel('$R_c$')
plt.ylabel('$N$')
plt.hist(slct_rcrit, bins_rc,
            histtype='step', color='g')
plt.tight_layout()

fontP = FontProperties()
fontP.set_size('xx-small')
fig_3   = plt.figure(figsize=(10,6))
ax_3    = fig_3.add_subplot(1, 1, 1)
for i in range(20):
    if nb_indexes[i] >= 10 :
   	plt.plot(med_profile[i,:], label="$R_{c} = $"+str(i))
plt.plot([0,150],[0.4,0.4])
plt.xlabel('$R [pix]$')
plt.ylabel('$Normalized Flux$')
plt.legend(loc=1, numpoints=1)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='small') 
plt.show()
