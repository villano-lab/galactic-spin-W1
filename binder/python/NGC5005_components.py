###############
### Imports ###
###############

import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm
import dataPython as dp
import scipy.integrate as si
from scipy.interpolate import InterpolatedUnivariateSpline

##############################
### Import data/text files ###
##############################

# Datapoints:
data = dp.getXYdata_wXYerr('data/NGC5005/ngc5005_data.txt')
r_dat = np.asarray(data['xx'])
v_dat = np.asarray(data['yy'])
v_err0 = np.asarray(data['ex'])
v_err1 = np.asarray(data['ey'])

# Disk:
disk_data = dp.getXYdata('data/NGC5005/ngc5005_disk.txt')
disk_raw = np.asarray(disk_data['yy'])

# Bulge:
bulge_data = dp.getXYdata('data/NGC5005/ngc5005_bulge.txt')
bulge_raw = np.asarray(bulge_data['yy'])

# Gas:
gas_data = dp.getXYdata('data/NGC5005/ngc5005_gas.txt')
gas_raw = np.asarray(gas_data['yy'])

#####################
### Interpolation ###
#####################

def interpd(x,y):
    return InterpolatedUnivariateSpline(x,y,k=5)

#######################
### Fitting Results ###
#######################

# Errors
weighdata = 1/v_err1

# Prefactors
bpref = 0.28228633
dpref = 0.71035496
gpref = 0.956

################################
######### Components ###########
################################

def bulge(r,bpref):
    polynomial = interpd(r_dat,bpref*bulge_raw)   
    return polynomial(r)

def disk(r,dpref):
    polynomial = interpd(r_dat,dpref*disk_raw)   
    return polynomial(r)

def gas(r,gpref):
    polynomial = interpd(r_dat,gpref*gas_raw)   
    return polynomial(r)


#########################
### Galaxy parameters ###
#########################

# NGC 5005 (Source: )
rc = 2.5                  # cutoff radius (in kpc) NEED REFERENCE
rho00 = 1.5981e+08         # central density (in solar mass/pc^3) NEED REFERENCE

# Constants
G = 4.30091e-6            # gravitational constant (kpc/solar mass*(km/s)^2)

#################################
### Calculating enclosed mass ### 
#################################

# NFW (dark halo) density profile
rho_NFW = lambda r: rho0 / ((r/rcut)*(1+r/rcut)**2)

# Inner function
mass_inner = lambda r: rho_NFW(r) * 4 * np.pi * r**2

# Mass integral: total mass at radius R (kpc)
#mass_r = lambda r: si.quad(mass_inner, 0, r)   
# Integral keeps giving me errors, so I used Mathematica to do the integral for me, this is the result:
def mass_r(r):
    return 4 * rcut**3 * np.pi * rho0 * (-1 + rcut/(rcut+r) - np.log(rcut) + np.log(rcut+r))

########################################################
### Calculating halo velocity using only black holes ###
########################################################

def halo_BH(r,mBH): 
    return mBH * np.sqrt(G * mass_r(r)/r)           # multiplied by a prefactor

##################################
### Calculating total velocity ###
##################################

def totalvelocity(r,mBH,bpref,dpref,gpref):
    return np.sqrt((disk(r,dpref))**2
               +(bulge(r,bpref))**2
               +(gas(r,gpref))**2
               +(halo_BH(r,mBH))**2)