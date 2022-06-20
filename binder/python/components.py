###############
### Imports ###
###############

import numpy as np
import dataPython as dp
import inspect
import types
import scipy.integrate as si
from scipy.interpolate import InterpolatedUnivariateSpline      # Spline function
import lmfit as lm                                              # Fitting
#custom libraries
import load_galaxies as galdata
import NGC5533_functions as funcNGC5533

ngc7814list = ['ngc7814', 'ngc 7814', '7814']
ngc5533list = ['ngc5533', 'ngc 5533', '5533']

##############################
### Import data/text files ###
##############################

# Note: there's no need to import the radius for each component as everything has the same r array 
# (the r array of the raw data)

# Datapoints:
def data(galaxy):
    if galaxy.lower() in ngc7814list:
        return dp.getXYdata_wYerr('data/NGC7814/7814_measured.dat')
    elif galaxy.lower() in ngc5533list:
        return galdata.data
def data_total(galaxy):
    if galaxy.lower() in ngc5533list:
        return galdata.data_total
    else:
        return None
def r_dat(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(data(galaxy)['xx'])
    elif galaxy.lower() in ngc5533list:
        return galdata.NGC5533['m_radii']
def v_dat(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(data(galaxy)['yy'])
    elif galaxy.lower() in ngc5533list:
        return galdata.NGC5533['m_velocities']
def v_err0(galaxy):
    if galaxy.lower() in ngc5533list:
        return galdata.v_err0
    else:
        return None
def v_err1(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(data(galaxy)['ey'])
    elif galaxy.lower() in ngc5533list:
        return galdata.NGC5533['m_v_errors']
    else:
        raise ValueError('Galaxy string', galaxy, 'not recognized.')
    
def disk_raw(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(dp.getXYdata('data/NGC7814/7814_gDisk.dat')['yy'])
    else:
        return None

def bulge_raw(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(dp.getXYdata('data/NGC7814/7814_gBulge.dat')['yy'])
    else:
        return None
    
def gas_raw(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(dp.getXYdata('data/NGC7814/7814_gGas.dat')['yy'])
    else:
        return None
    
#####################
### Interpolation ###
#####################

def interpd(x,y):
    return InterpolatedUnivariateSpline(x,y,k=5)

################################
######### Components ###########
################################

def blackhole(r,M,galaxy):
    if galaxy.lower() in ngc5533list:
        x = np.sort(r_dat(galaxy))
        y = funcNGC5533.bh_v(r_dat(galaxy),M,load=False)
        polynomial = interpd(x,y)
        return polynomial(r)
    else: #galaxy has no significant bh component
        return 0
    
def bulge(r,bpref,galaxy):
    if galaxy.lower() in ngc7814list:
        polynomial = interpd(r_dat(galaxy),bpref*bulge_raw(galaxy))   
        return polynomial(r)
    elif galaxy.lower() in ngc5533list:
        x = np.sort(r_dat(galaxy))
        y = bpref*funcNGC5533.b_v(r_dat(galaxy),load=True)
        polynomial = interpd(x,y)
        return polynomial(r)

def disk(r,dpref,galaxy):
    if galaxy.lower() in ngc7814list:
        polynomial = interpd(r_dat(galaxy),dpref*disk_raw(galaxy))   
        return polynomial(r)
    elif galaxy.lower() in ngc5533list:
        return funcNGC5533.d_thief(r,dpref)

def gas(r,gpref,galaxy):
    if galaxy.lower() in ngc7814list:
        polynomial = interpd(r_dat(galaxy),gpref*gas_raw(galaxy))   
        return polynomial(r)
    elif galaxy.lower() in ngc5533list:
        return funcNGC5533.g_thief(r,gpref)
    
#########################
### Galaxy parameters ###
#########################


def rcut(galaxy):
    if galaxy.lower() in ngc7814list: # NGC 7814 (Source: Fraternali, Sancisi, and Kamphuis, 2011)
        return 2.1            # Cutoff (core) radius (in kpc) from Table 5.
    elif galaxy.lower() in ngc5533list:
        return funcNGC5533.h_rc
    
def rho0(galaxy):
    if galaxy.lower() in ngc7814list:
        return 1.524e8            # Central density (in solar mass/kpc^3) from Table 5. (converted from pc to kpc)
    elif galaxy.lower() in ngc5533list:
        return funcNGC5533.hrho00_c

# Constants
G = 4.30091e-6            # Gravitational constant (kpc/solar mass*(km/s)^2)

#################################
### Calculating enclosed mass ### 
#################################

# Mass as a function of radius, calculated from the Isothermal density profile:
def mass_r(r,rcut):
    return 4 * np.pi * rcut**2  * (r - rcut * (np.arctan(r/rcut))) 
# rho0 represents the number of tiny black holes at the center of the galaxy

###############################################################
### Calculating halo velocity normally or using black holes ###
###############################################################

# Calculating the velocity for each black hole as a point mass
def halo_BH(r,scale,arraysize,massMiniBH,rcut):
    x = np.sort(r)
    y = np.sqrt((G * (scale * arraysize * massMiniBH) * mass_r(r,rcut)) / r)   
                    # scale is needed to be separate and constant because the widget would freeze the computer otherwise
                    # arraysize is the number of black holes for slider
                    # mBH is the mass of black holes for slider
    halo = interpd(x,y)
    return halo(r)

def halo(r,rc,rho00,galaxy):
    if galaxy in ngc5533list:
        return funcNGC5533.h_v(r,rc,rho00,load=False)

##################################
### Calculating total velocity ###
##################################
def totalvelocity_miniBH(r,scale,arraysize,massMiniBH,rcut,bpref,dpref,gpref,Mbh,galaxy):
    return np.sqrt(blackhole(r,Mbh,galaxy)**2 
                        + bulge(r,bpref,galaxy)**2 
                        + disk(r,dpref,galaxy)**2
                        + halo_BH(r,scale,arraysize,massMiniBH,rcut)**2
                        + gas(r,gpref,galaxy)**2)
    
def totalvelocity_halo(r,scale,arraysize,rho00,rcut,bpref,dpref,gpref,Mbh,galaxy):
    return np.sqrt(funcNGC5533.bh_v(r,Mbh,load=False)**2
                   + bulge(r,bpref,galaxy)**2
                   + disk(r,dpref,galaxy)**2
                   + halo(r,rcut,rho00,galaxy)**2
                   + gas(r,gpref,galaxy)**2)
    
#################################
#### Find Fitting Parameters ####
#################################
# For tiny black hole widget

# Setup
def weighdata(galaxy):
    return 1/v_err1(galaxy)

# Tiny Black Holes 
def set_params(model,galaxy):
    if type(model) == types.FunctionType:
        model = lm.Model(model)
        fit_pars = model.make_params()
    elif type(model) == lm.Model:
        fit_pars = model.make_params()
    else:
        raise ValueError("Invalid type for variable `model`. (",model,", type:",type(model),".)")
    #miniBH halo
    #fit_pars.add('rc',    value=rcut(galaxy), min=0.1)         # Core Radius (kpc)
    fit_pars.add('scale',      value=rho0(galaxy),   vary=False)        # Scale
    if (model == 'bh') or (model==lm.Model(totalvelocity_miniBH)) or (model==totalvelocity_miniBH):
        fit_pars.add('arraysize',  value=50,     min=1, max=100)    # Number of black holes
        fit_pars.add('rho0', value=1.5, min=0)       # Halo Central Density 
    elif (model == 'wimp') or (model==lm.Model(totalvelocity_halo)) or (model==totalvelocity_halo):
        fit_pars.add('arraysize', value=0, vary=False)
        fit_pars.add('rho0', value=funcNGC5533.hrho00_c, min=0)
    fit_pars.add('rcut',       value=rcut(galaxy),   min=0.1)           # Core Radius (kpc)
    # Bulge
    fit_pars.add('bpref', value=1, min=0, max=100)  # Bulge Prefactor
    # Disk
    fit_pars.add('dpref', value=1, min=0, max=100)  # Disk Prefactor
    # Disk
    fit_pars.add('gpref', value=1, vary=False)        # Gas Prefactor
    # BH
    if galaxy in ngc5533list:
        fit_pars.add('Mbh', value=funcNGC5533.Mbh_def, min=1e8)
    else:
        fit_pars.add('Mbh', value=0, vary=False)
    return fit_pars

# Do fit
def bestfit(model,galaxy):
    newmodel = lambda r,scale,arraysize,rho0,rcut,bpref,dpref,gpref,Mbh: model(r,scale,arraysize,rho0,rcut,bpref,dpref,gpref,Mbh,galaxy)
    fit_mod = lm.Model(newmodel)
    fit_pars = set_params(newmodel,galaxy)
    if (model == 'bh') or (model==lm.Model(totalvelocity_miniBH)) or (model==totalvelocity_miniBH):
        fit_pars.add('arraysize',  value=50,     min=1, max=100)    # Number of black holes
        fit_pars.add('rho0', value=1.5, min=0)       # Halo Central Density 
    elif (model == 'wimp') or (model==lm.Model(totalvelocity_halo)) or (model==totalvelocity_halo):
        fit_pars.add('arraysize', value=0, vary=False)
        fit_pars.add('rho0', value=funcNGC5533.hrho00_c, min=0)
    if galaxy.lower() in ngc5533list and (model == lm.Model(totalvelocity_halo) or model == totalvelocity_halo):
        weights = 1/np.sqrt(v_err1(galaxy)**2+galdata.NGC5533['n_v_bandwidth']**2)
    else:
        weights = weighdata(galaxy)
    fit = fit_mod.fit(v_dat(galaxy),fit_pars,r=r_dat(galaxy),weights=weights)
    bestfit = fit.best_fit
    fit_dict = fit.best_values
    return bestfit, fit_dict
    #Examples on how to use this output:
    #best_rc = fit_dict['rc']
    #best_rho00 = fit_dict['rho00']
    #best_bpref = fit_dict['bpref']
    #best_dpref = fit_dict['dpref']
    #best_gpref = fit_dict['gpref']