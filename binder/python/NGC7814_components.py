###############
### Imports ###
###############

import numpy as np
import dataPython as dp
import scipy.integrate as si
from scipy.interpolate import InterpolatedUnivariateSpline      # Spline function
import lmfit as lm                                              # Fitting

##############################
### Import data/text files ###
##############################

# Note: there's no need to import the radius for each component as everything has the same r array 
# (the r array of the raw data)

# Datapoints:
data = dp.getXYdata_wXYerr('data/NGC7814/7814_measured.dat')
r_dat = np.asarray(data['xx'])
v_dat = np.asarray(data['yy'])
v_err1 = np.asarray(data['ey'])

# Disk:
disk_data = dp.getXYZdata('data/NGC7814/7814_gDisk.dat')
disk_raw = np.asarray(disk_data['zz'])

# Bulge:
bulge_data = dp.getXYZdata('data/NGC7814/7814_gBulge.dat')
bulge_raw = np.asarray(bulge_data['zz'])

# Gas:
gas_data = dp.getXYZdata('data/NGC7814/7814_gGas.dat')
gas_raw = np.asarray(gas_data['zz'])

#####################
### Interpolation ###
#####################

def interpd(x,y):
    return InterpolatedUnivariateSpline(x,y,k=5)

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

# NGC 7814 (Source: Fraternali, Sancisi, and Kamphuis, 2011)
rcut = 2.1                # Cutoff (core) radius (in kpc) from Table 5.
rho0 = 1.524e8            # Central density (in solar mass/kpc^3) from Table 5. (converted from pc to kpc)

# Constants
G = 4.30091e-6            # Gravitational constant (kpc/solar mass*(km/s)^2)

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
def halo_BH(r,scale,arraysize,mBH,rcut):
    x = np.sort(r)
    y = np.sqrt((G * (scale * arraysize * mBH) * mass_r(r,rcut)) / r)   
                    # scale is needed to be separate and constant because the widget would freeze the computer otherwise
                    # arraysize is the number of black holes for slider
                    # mBH is the mass of black holes for slider
    halo = interpd(x,y)
    return halo(r)

##################################
### Calculating total velocity ###
##################################

def totalvelocity(r,scale,arraysize,massMiniBH,rcut,bpref,dpref,gpref):    
    total = np.sqrt(bulge(r,bpref)**2 
                    + disk(r,dpref)**2
                    + halo_BH(r,scale,arraysize,massMiniBH,rcut)**2
                    + gas(r,gpref)**2)
    return total

#################################
#### Find Fitting Parameters ####
#################################
# For tiny black hole widget

# Setup
fit_mod = lm.Model(totalvelocity)
fit_params = fit_mod.make_params()
weighdata = 1/v_err1

# Tiny Black Holes 
fit_params.add('scale',      value=2.5e7,  vary=False)        # Scale
fit_params.add('arraysize',  value=50,     min=1, max=100)    # Number of black holes
fit_params.add('massMiniBH', value=1.5,    min=0.1, max=3.8)  # Mass of each tiny black holes (Solar mass)
fit_params.add('rcut',       value=rcut,   min=0.1)           # Core Radius (kpc)
# Bulge
fit_params.add('bpref', value=1, min=0.5, max=100)  # Bulge Prefactor
# Disk
fit_params.add('dpref', value=1, min=0.5, max=100)  # Disk Prefactor
# Disk
fit_params.add('gpref', value=1, vary=False)        # Gas Prefactor

# Do fit
fit = fit_mod.fit(v_dat,fit_params,r=r_dat,weights=weighdata)

#################################
### Define fitting parameters ###
#################################

bestfit = fit.best_fit

fit_dict = fit.best_values
best_scale = fit_dict['scale']
best_arraysize = fit_dict['arraysize']
best_massMiniBH = fit_dict['massMiniBH']
best_rcut = fit_dict['rcut']
best_bpref = fit_dict['bpref']
best_dpref = fit_dict['dpref']
best_gpref = fit_dict['gpref']