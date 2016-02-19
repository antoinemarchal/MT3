import matplotlib.pyplot as plt
import numpy as np
import pickle

slct_redshift   = pickle.load(open("results/slct_redshift.pkl","rb"))
slct_flux       = pickle.load(open("results/slct_flux.pkl","rb"))
slct_msz        = pickle.load(open("results/slct_msz.pkl","rb"))
RD_in_redshift  = pickle.load(open("results/RD_in_redshift.pkl","rb")) 
RD_out_redshift = pickle.load(open("results/RD_out_redshift.pkl","rb"))
RD_in_flux      = pickle.load(open("results/RD_in_flux.pkl","rb"))
RD_out_flux     = pickle.load(open("results/RD_out_flux.pkl","rb"))
RD_in_msz       = pickle.load(open("results/RD_in_msz.pkl","rb"))
RD_out_msz      = pickle.load(open("results/RD_out_msz.pkl","rb"))
RD_in_rcrit     = pickle.load(open("results/RD_in_rcrit.pkl","rb"))
RD_out_rcrit    = pickle.load(open("results/RD_out_rcrit.pkl","rb"))
slct_rcrit      = pickle.load(open("results/slct_rcrit.pkl","rb"))

moy   = np.mean(slct_flux)
std   = np.std(slct_flux)
moy_msz=np.mean(RD_in_msz)
"""----------Plot/Results Study flow/redshift------------"""
bins_r = np.linspace(0., np.max(RD_in_redshift), 40)
bins_f = np.linspace(np.min(RD_in_flux),
                     np.max(RD_in_flux), 40)
bins_m = np.linspace(np.min(slct_msz),
                     np.max(slct_msz), 40)
bins_rc = np.linspace(np.min(slct_rcrit),
                     np.max(slct_rcrit), 40)
color = ['b.','g.', 'r.', 'c.', 'm.', 'k.', 'y.'] 
fig_1   = plt.figure()#figsize=(16,9))
ax_1    = fig_1.add_subplot(1, 1, 1)
plt.locator_params(nbins=4)
plt.subplot(2,2,1)
plt.xlabel('$z$')
plt.ylabel('$Flux$')
plt.plot(RD_in_redshift, RD_in_flux, color[2],
         RD_out_redshift, RD_out_flux, color[6])

plt.subplot(2,2,2)
plt.xlabel('$z$')
plt.ylabel('$Flux$')
plt.plot(RD_in_redshift, RD_in_flux, color[2])

plt.subplot(2,2,3)
plt.xlabel('$Flux$')
plt.ylabel('$N_{cluster}$')
plt.hist(RD_in_flux, bins_f,
            histtype='step', color='g')
plt.plot([moy, moy], [0,100], 'r--', lw=2)

plt.subplot(2,2,4)
plt.xlabel('$z$')
plt.ylabel('$N_{cluster}$')
plt.hist(RD_in_redshift, bins_r,
            histtype='step', color='b')
plt.tight_layout()
plt.savefig('results/flux_z.pdf', format='pdf')

fig_2   = plt.figure()#figsize=(16,9))
ax_2    = fig_2.add_subplot(1, 1, 1)
plt.subplot(2,2,1)
plt.xlabel('$M/10^{14} M \odot$')
plt.ylabel('$Flux$')
plt.plot(RD_in_msz, RD_in_flux, color[4],
         RD_out_msz, RD_out_flux, color[6])

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
plt.xlabel('$r_c$')
plt.ylabel('$N$')
plt.hist(slct_rcrit, bins_rc,
            histtype='step', color='g')
plt.tight_layout()

plt.show()
