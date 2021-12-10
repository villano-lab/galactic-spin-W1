###############
### Imports ###
###############

import numpy as np
import dataPython as dp
import scipy.integrate as si
from scipy.interpolate import InterpolatedUnivariateSpline

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

#######################
### Fitting Results ###
#######################

# Errors
weighdata = 1/v_err1

# LMFit results from fitting
best_bpref = 4.98331928       # bulge
best_dpref = 3.85244293       # disk
best_gpref = 1.00239573       # gas

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

# NGC 7814 (Source: https://www.aanda.org/articles/aa/abs/2011/07/aa16634-11/aa16634-11.html)
rcut = 2.1                # cutoff radius (in kpc) from Table 5.
rho0 = 152.3e-3           # central density (in solar mass/pc^3) from Table 5.

# Constants
G = 4.30091e-6            # gravitational constant (kpc/solar mass*(km/s)^2)

#################################
### Calculating enclosed mass ### 
#################################

# Mass as a function of radius, calculated from the Isothermal density profile (see Mathematica code):
def mass_r(r,rcut):
    return 4 * np.pi * rcut**2  * (r - rcut * (np.arctan(r/rcut)))    

########################################################
### Calculating halo velocity using only black holes ###
########################################################

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

def totalvelocity(r,scale,arraysize,mBH,rcut,bpref,dpref,gpref):    
    total = np.sqrt(bulge(r,bpref)**2 
                    + disk(r,dpref)**2
                    + halo_BH(r,scale,arraysize,mBH,rcut)**2
                    + gas(r,gpref)**2)
    return total