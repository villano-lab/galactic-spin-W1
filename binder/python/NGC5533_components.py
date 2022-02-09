###############
### Imports ###
###############

import numpy as np
import dataPython as dp
import scipy.integrate as si
from scipy.interpolate import InterpolatedUnivariateSpline

import NGC5533_functions as funcNGC5533                   # Functions for the components of NGC 5533
import NGC5533_traced_curves as dataNGC5533               # Traced data of Galaxy NGC 5533
import NGC5533_fitting as fitNGC5533                      # Fitting parameters for best fit values

################################
############ Data ##############
################################

data = dataNGC5533.data
data_total = dataNGC5533.data_total
r_dat = dataNGC5533.r_dat
v_dat = dataNGC5533.v_dat
v_err0 = dataNGC5533.v_err0
v_err1 = dataNGC5533.v_err1

#####################
### Interpolation ###
#####################

def interpd(x,y):
    return InterpolatedUnivariateSpline(x,y,k=5)

################################
######### Components ###########
################################

def blackhole(r,M):
    x = np.sort(r_dat)
    y = funcNGC5533.bh_v(r_dat,M,load=False)
    polynomial = interpd(x,y)
    return polynomial(r)

def bulge(r,bpref):
    x = np.sort(r_dat)
    y = bpref*funcNGC5533.b_v(r_dat,load=True)
    polynomial = interpd(x,y)
    return polynomial(r)

def disk(r,dpref):
    return funcNGC5533.d_thief(r,dpref)

def gas(r,gpref):
    return funcNGC5533.g_thief(r,gpref)

################################
###### Fitting Parameters ######
################################

best_Mbh = fitNGC5533.f_Mbh
best_bpref = fitNGC5533.f_bpref
best_dpref = fitNGC5533.f_dpref
#rcut = fitNGC5533.f_rc      # value determined by slider
#rho0 = fitNGC5533.f_rho00  # don't need rho0, we're dealing with mass
best_gpref = fitNGC5533.f_gpref

# Constants
G = 4.30091e-6            # gravitational constant (kpc/solar mass*(km/s)^2)

#################################
### Calculating enclosed mass ### 
#################################

# Mass as a function of radius, calculated from the Isothermal density profile:
def mass_r(r,rcut):
    return 4 * np.pi * rcut**2  * (r - rcut * (np.arctan(r/rcut)))    
# rho0 represents the number of tiny black holes at the center of the galaxy

########################################################
### Calculating halo velocity using only black holes ###
########################################################

# Calculating the velocity for each black hole as a point mass
def halo_BH(r,scale,arraysize,massMiniBH,rcut):
    x = np.sort(r)
    y = np.sqrt((G * (scale * arraysize * massMiniBH) * mass_r(r,rcut)) / r)   
                    # scale is needed to be separate and constant because the widget would freeze the computer otherwise
                    # arraysize is the number of black holes for slider
                    # mBH is the mass of black holes for slider
    halo = interpd(x,y)
    return halo(r)

##################################
### Calculating total velocity ###
##################################

def totalvelocity(r,scale,arraysize,massMiniBH,rcut,Mbh,bpref,dpref,gpref):    
    total = np.sqrt(blackhole(r,Mbh)**2 
                    + bulge(r,bpref)**2 
                    + disk(r,dpref)**2
                    + halo_BH(r,scale,arraysize,massMiniBH,rcut)**2
                    + gas(r,gpref)**2)
    return total