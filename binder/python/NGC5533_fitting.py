################################
########### Imports ############
################################

import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm
import dataPython as dp
import scipy.interpolate as inter

import NGC5533_functions as funcNGC5533                   # Functions for the components of NGC 5533
import NGC5533_traced_curves as dataNGC5533               # Traced data of Galaxy NGC 5533
import NGC5533_fitting as fitNGC5533                      # Fitting parameters for best fit values

################################
##### Measured data points #####
################################

data = dp.getXYdata_wXYerr('data/NGC5533/100kpc_data.txt')
r_dat = np.asarray(data['xx'])
v_dat = np.asarray(data['yy'])
v_err0 = np.asarray(data['ex'])
v_err1 = np.asarray(data['ey'])

###########################################
###### Uncertainty Band and Weights #######
###########################################

# Uncertainty band, using Noordermeer's function with a guessed delta_i
delta_i = 3     # guessed value
v_i = (v_dat / np.tan(52*(np.pi/180)) * delta_i *(np.pi/180))

# Traced uncertainty band
band = dataNGC5533.band
greyb_bottom = dataNGC5533.greyb_bottom
greyb_top = dataNGC5533.greyb_top

# Change r_dat so it's strictly increasing
r_dat, v_dat, v_err0, v_err1 = (np.asarray(list(a)) for a in zip(*sorted(zip(r_dat, v_dat, v_err0, v_err1))))
# Converting v_err1 to an array
v_err1_array = np.asarray(v_err1)
# Express as weights
#weighdata = 1/(np.sqrt((v_err1**2)+(v_i**2)))
weighdata = 1/(np.sqrt((v_err1**2)+(band**2)))

################################
########## Function ############
################################

r = np.arange(0.1,200,0.1)

def f(r,Mbh,rc,rho00,bpref,dpref,gpref):
    return np.sqrt(funcNGC5533.bh_v(r,Mbh,load=False)**2 
                   + funcNGC5533.h_v(r,rc,rho00,load=False)**2 
                   + bpref**2*funcNGC5533.b_v(r,load=True,path='../')**2 
                   + dpref**2*funcNGC5533.d_thief(r)**2
                   + gpref**2*funcNGC5533.g_thief(r)**2)

##########################################
########## Fitting Parameters ############
##########################################

# Setup
f_mod = lm.Model(f)
f_params = f_mod.make_params()
# Black Hole
f_params.add('Mbh',   value=funcNGC5533.Mbh_def, min=1.0e8)    # Black Hole Mass
# Halo
f_params.add('rc',    value=funcNGC5533.h_rc, min=0.1)         # Core Radius (kpc)
f_params.add('rho00', value=funcNGC5533.hrho00_c, min=0)       # Halo Central Density 
# Bulge
f_params.add('bpref', value=1,min=0,max=100)                   # Bulge Prefactor
# Disk
f_params.add('dpref', value=1,min=0, max=100)                  # Disk Prefactor
# Gas
f_params.add('gpref', value=1, vary=False)                     # Gas Prefactor

# Do fit
f_fit = f_mod.fit(v_dat,f_params,r=r_dat,weights=weighdata)

##########################################
########## Define for Plotting ###########
##########################################

bestf = f_fit.best_fit
#delf = f_fit.eval_uncertainty()

f_dict = f_fit.best_values
f_Mbh = f_dict['Mbh']
f_rc = f_dict['rc']
f_rho00 = f_dict['rho00']
f_bpref = f_dict['bpref']
f_dpref = f_dict['dpref']
f_gpref = f_dict['gpref']

f_curve = f(r,f_Mbh,f_rc,f_rho00,f_bpref,f_dpref,f_gpref)
bh_curve = funcNGC5533.bh_v(r,f_Mbh,load=False)
halo_curve = funcNGC5533.h_v(r,f_rc,f_rho00,load=False)
bulge_curve = f_bpref*funcNGC5533.b_v(r,load=True)
disk_curve = f_dpref*funcNGC5533.d_thief(r)
gas_curve = f_gpref*funcNGC5533.g_thief(r)
