################################
########### Imports ############
################################

import sys 
#sys.path.append('../python/')
import NGC5533_functions as nf
import noordermeer as noord               # Traced data of Galaxy NGC 5533
import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm
import dataPython as dp
import scipy.interpolate as inter

################################
##### Measured data points #####
################################

data = dp.getXYdata_wXYerr('../data/100kpc_data.txt')
#data = dp.getXYdata_wXYerr('../galactic-spin-binder/NGC_5533/data/final/100kpc_data.txt')
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
band = noord.band
greyb_bottom = noord.greyb_bottom
greyb_top = noord.greyb_top

#change r_dat so it's strictly increasing
r_dat, v_dat, v_err0, v_err1 = (np.asarray(list(a)) for a in zip(*sorted(zip(r_dat, v_dat, v_err0, v_err1))))
#converting v_err1 to an array
v_err1_array = np.asarray(v_err1)
# Express as weights
#weighdata = 1/(np.sqrt((v_err1**2)+(v_i**2)))
weighdata = 1/(np.sqrt((v_err1**2)+(band**2)))

################################
########## Function ############
################################

r = np.arange(0.1,200,0.1)

def f(r,M,rc,rho00,c,pref,gpref):
    return np.sqrt(nf.bh_v(r,M,load=False)**2 
                   + nf.h_v(r,rc,rho00,load=False)**2 
                   + c**2*nf.b_v(r,load=True,path='../')**2 
                   + pref**2*nf.d_thief(r)**2
                   + gpref**2*nf.g_thief(r)**2)

##########################################
########## Fitting Parameters ############
##########################################

# Setup
f_mod = lm.Model(f)
f_params = f_mod.make_params()
# Black Hole
f_params.add('M', value=nf.Mbh_def, min=1.0e8)     #Mass
# Halo
f_params.add('rc', value=nf.h_rc, min=0.1)         #Core Radius (kpc)
f_params.add('rho00', value=nf.hrho00_c, min=0)    #Halo Density 
# Bulge
f_params.add('c', value=1,min=0,max=100)           #Prefactor
# Disk
f_params.add('pref', value=1,min=0, max=100)       #Prefactor
# Gas
f_params.add('gpref', value=1, vary=False)     #Prefactor

# Do fit
f_fit = f_mod.fit(v_dat,f_params,r=r_dat,weights=weighdata)

##########################################
########## Define for Plotting ###########
##########################################

bestf = f_fit.best_fit
#delf = f_fit.eval_uncertainty()

f_dict = f_fit.best_values
f_M = f_dict['M']
f_c = f_dict['c']
f_pref = f_dict['pref']
f_rc = f_dict['rc']
f_hrho00 = f_dict['rho00']
f_gpref = f_dict['gpref']

f_curve = f(r,f_M,f_rc,f_hrho00,f_c,f_pref,f_gpref)
bh_curve = nf.bh_v(r,f_M,load=False)
halo_curve = nf.h_v(r,f_rc,f_hrho00,load=False)
bulge_curve = f_c*nf.b_v(r,load=True)
disk_curve = f_pref*nf.d_thief(r)
gas_curve = f_gpref*nf.g_thief(r)
