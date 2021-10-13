################################
########### Imports ############
################################

import sys 
#sys.path.append('../python/')
import NGC5533_functions as nf
import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm
import dataPython as dp
import scipy.interpolate as inter
from scipy.interpolate import InterpolatedUnivariateSpline

#######################
##### Traced data #####
#######################

#data
data = dp.getXYdata_wXYerr('testing/891/891_data')
r_dat = np.asarray(data['xx'])
v_dat = np.asarray(data['yy'])
v_err1 = np.asarray(data['ey'])

#disk:
disk = dp.getXYdata('testing/891/891_dtDisk.dat')
r_d = np.asarray(disk['xx'])
v_d = np.asarray(disk['yy'])

#bulge:
bulge = dp.getXYdata('testing/891/891_dtBulge.dat')
r_b = np.asarray(bulge['xx'])
v_b = np.asarray(bulge['yy'])

#gas:
gas = dp.getXYdata('testing/891/891_dtGas.dat')
r_g = np.asarray(gas['xx'])
v_g = np.asarray(gas['yy'])

#####################
### Interpolation ###
#####################

def interpd(x,y):
    return InterpolatedUnivariateSpline(x,y,k=4)

################################
######### Components ###########
################################

def bulge(r):
    x = r_b
    y = v_b
    polynomial = interpd(x,y)
    return polynomial(r)

def disk(r):
    x = r_d
    y = v_d
    polynomial = interpd(x,y)
    return polynomial(r)

def gas(r):
    x = r_g
    y = v_g
    polynomial = interpd(x,y)
    return polynomial(r)

######################
###### Weights #######
######################

v_err1 = v_err1
weighdata = 1/v_err1

##########################
### Default parameters ###
##########################

rcut = 2.5      # cutoff radius (in kpc)
rho0 = 1e8      # central density

################################
########## Function ############
################################

r = np.linspace(0.1,16.4,1000)

def t(r,bpref,dpref,rcut,rho0,gpref):
    return np.sqrt((bpref*bulge(r))**2
                   +(dpref*disk(r))**2
                   +(nf.h_v(r,rcut,rho0))**2
                   +(gpref*gas(r))**2)

##########################################
########## Fitting Parameters ############
##########################################

# Setup
g_mod = lm.Model(t)
g_params = g_mod.make_params()
# Halo
g_params.add('rcut', value=rcut, min=0.1)          #Core radius (kpc)
g_params.add('rho0', value=rho0, min=0)            #Central density 
# Bulge
g_params.add('bpref', value=1, min=0)              #Prefactor
# Disk
g_params.add('dpref', value=1, min=0)              #Prefactor
# Gas
g_params.add('gpref', value=1, min=.99,max=1.01)   #Prefactor
# Do fit
g_fit = g_mod.fit(v_dat,g_params,r=r_dat,weights=weighdata)

##########################################
########## Define for Plotting ###########
##########################################

bestf = g_fit.best_fit

g_dict = g_fit.best_values
best_bpref = g_dict['bpref']
best_dpref = g_dict['dpref']
best_rcut = g_dict['rcut']
best_rho0 = g_dict['rho0']
best_gpref = g_dict['gpref']