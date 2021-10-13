###############
### Imports ###
###############

import numpy as np
import sys
sys.path.append('../../python/')
import dataPython as dp
import scipy.integrate as si
from scipy.interpolate import InterpolatedUnivariateSpline

import NGC5533_functions as nf            # Components
import noordermeer as noord               # Traced data of Galaxy NGC 5533
import fitting_NGC5533 as fitting         # Fitting parameters for best fit values

################################
############ Data ##############
################################

data = noord.data
data_total = noord.data_total
r_dat = noord.r_dat
v_dat = noord.v_dat
v_err0 = noord.v_err0
v_err1 = noord.v_err1

#####################
### Interpolation ###
#####################

def interpd(x,y):
    return InterpolatedUnivariateSpline(x,y,k=5)

################################
######### Components ###########
################################

def blackhole(r,M):
    x = np.sort(r)
    y = nf.bh_v(r,M,load=False)
    polynomial = interpd(x,y)
    return polynomial(r)

def bulge(r,bpref):
    x = np.sort(r)
    y = bpref*nf.b_v(r,load=True)
    polynomial = interpd(x,y)
    return polynomial(r)

def disk(r,dpref):
    x = np.sort(r)
    y = dpref*nf.d_thief(r)
    polynomial = interpd(x,y)
    return polynomial(r)

def gas(r,gpref):
    x = np.sort(r)
    y = gpref*nf.g_thief(r)
    polynomial = interpd(x,y)
    return polynomial(r)

################################
###### Fitting Parameters ######
################################

best_M = fitting.f_M
best_bpref = fitting.f_c
best_dpref = fitting.f_pref
#rcut = fitting.f_rc      # value determined by slider
#rho0 = fitting.f_hrho00  # don't need rho0, we're dealing with mass
best_gpref = fitting.f_gpref

# Constants
G = 4.30091e-6            # gravitational constant (kpc/solar mass*(km/s)^2)

#################################
### Calculating enclosed mass ### 
#################################

# Mass as a function of radius, calculated from the Isothermal density profile (see Mathematica code):
def mass_r(r,rcut):
    return 4 * np.pi * rcut**2  * (r - rcut * (np.arctan(r/rcut)))    
# rho0 represents the number of tiny black holes at the center of the galaxy

########################################################
### Calculating halo velocity using only black holes ###
########################################################

# What if I just calculate the velocity for each black hole as a point mass?
def halo_BH(r,scale,arraysize,mBH,rcut):
    x = np.sort(r)
    y = np.sqrt((G * (scale * arraysize * mBH) * mass_r(r,rcut)) / r)   
                    # scale is needed to be separate and constant because the widget would freeze the computer otherwise
                    # arraysize is the number of black holes for slider
                    # mBH is the mass of black holes for slider
    polynomial = interpd(x,y)
    return polynomial(r)

##################################
### Calculating total velocity ###
##################################

def totalvelocity(r,scale,arraysize,mBH,rcut,M,bpref,dpref,gpref):    
    total = np.sqrt(blackhole(r,M)**2 
                    + bulge(r,bpref)**2 
                    + disk(r,dpref)**2
                    + halo_BH(r,scale,arraysize,mBH,rcut)**2
                    + gas(r,gpref)**2)
    return total