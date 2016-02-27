#SZ_effect
Data analysis / Planck 2015
Using ILC method to extract thermal Sunyaev–Zel'dovich effect maps from Planck HFI data.

Use HFI maps (2048) of Planck mission and the PSZ2v1.fits catalog available on SZ Cluster DataBase (IAS/IDOC).

HFI_SkyMap_*.fits --> /SZ_effect/maps_2015
PSZ2v1.fits  --> /SZ_effect

--> Run smooth_psf.py 
    to smooth each map as a function of the worst resolution (100 GHz).

--> Run find_SZE.py [-h] [-p] n_cluster
    to extract patch map from the data where the PSZ2 catalog indicate the presence of galaxy cluster.
    We select the clusters based on the weight/latidude function and we save *.fits in patch_SZ/SZ directory.  
    (see report for details)

--> Run apert_phot.py 
    to applied aperture photometry on the extracted patchs and compute the radial profiles and flux of each cluster.
    The different plots are saved in the /results directory.
    
    
