from astropy.io import fits as pyfits
import pyfits as pf
from astropy import wcs
import astropy.table as pytabs
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import math as ma
import os
import function as fct
import mod_ilc as ilc
plt.ion()

def dist_SZ(nu):
    """----------------------------------------------------
    --- fact_SZE : Getting back the factor f(nu) 
    
                   Parameters : nu  = frequency in GHz
    
                   Return     : f   = adimenssionless factor
    --------------------------------------------------------"""
    import astropy.units as u
    from astropy import constants as const
    x = const.h * nu * u.GHz / const.k_B / ( 2.725 * u.K)
    f = x * ((np.exp(x) + 1.) / (np.exp(x) - 1.)) - 4.
    return f

def weight(R, a, a_t, b, b_t):
    """-------------------------------------------------------
    --- weight : Compute the weight of each map to obtein
                 CMB+KSZ map and TSZ map.
    
                 Parameters : R      = matrix of covariance
                                       between the maps
                              a, a_t = (1,...,1) and T matrix
                              b, b_t = f(nu) and T matrix
    
                 Return     : w_k    = weight vector of CMB+KSZ 
                                       map
                              w_t    = weight vector of TSZ map
    -----------------------------------------------------------"""
    g_1 = np.matmul(np.matmul(b_t,R),b) * (np.matmul(a_t,R))
    g_2 = np.dot(np.matmul(a_t,R),b) * (np.matmul(b_t,R))
    g_3 = np.matmul(np.matmul(a_t,R),b) * np.matmul(np.matmul(b_t,R),b)
    g_4 = np.matmul(np.matmul(a_t,R),b)**2
    
    f_1 = np.matmul(np.matmul(a_t,R),a) * (np.matmul(b_t,R))
    f_2 = np.matmul(np.matmul(b_t,R),a) * (np.matmul(a_t,R))
    f_3 = np.matmul(np.matmul(b_t,R),b) * np.matmul(np.matmul(a_t,R),a)
    f_4 = np.matmul(np.matmul(b_t,R),a)**2
    
    w_k   = (g_1 - g_2) / (g_3 - g_4)
    w_t   = (f_1 - f_2) / (f_3 - f_4)
    return (w_k, w_t)

def separation(k, map_smooth, patch_size, NAME, GLON, GLAT, freq):
    """------------------------------------------------------------
    --- separation :  
    
                     Parameters : FIXME
    
                     Return     : w_k    = weight vector of CMB+KSZ 
                                           map
                                  w_t    = weight vector of TSZ map
    ---------------------------------------------------------------"""
    n_obs         = 6
    f_nu          = ilc.dist_SZ(freq)
    
    a             = np.ones(n_obs)
    a_t           = np.transpose(a)
    b             = np.ones(n_obs)
    b             = f_nu
    b_t           = np.transpose(b)
    
    patch_map     = []
    
    for i in range(6):
        patch_map.append(
            (fct.patch_map(map_smooth[i], patch_size, GLON[k], GLAT[k]))
        ) 
        
    E = np.cov((patch_map[0].flatten(), patch_map[1].flatten(),
                patch_map[2].flatten(), patch_map[3].flatten(),
                patch_map[4].flatten(), patch_map[5].flatten()))
    R = np.linalg.inv(E)
        
    w_k, w_t = ilc.weight(R, a, a_t, b ,b_t)
        
    CMB_KSZ_map   = np.zeros((patch_size, patch_size))
    TSZ_map       = np.zeros((patch_size, patch_size))
        
    for i in range(n_obs) :
        CMB_KSZ_map = CMB_KSZ_map + (w_k[i]    * patch_map[i])
        TSZ_map     = TSZ_map     + (w_t[i]    * patch_map[i])
    return (CMB_KSZ_map, TSZ_map, w_t, w_k)

    
