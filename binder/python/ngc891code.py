#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm
import sys
sys.path.append('../python')
import dataPython as dp
import NGC5533_functions as nf


# In[2]:


#**********************importing text files******************************
#there's no need to import the radius for each component as everything has the same r array (the r array of the raw data)
#data:
data = dp.getXYdata_wXYerr('../testing/891/891_data')
r = np.asarray(data['xx'])
v_dat = np.asarray(data['yy'])
v_err1 = np.asarray(data['ey'])
#disk:
disk = dp.getXYdata('../testing/891/891_dtDisk.dat')
d = np.asarray(disk['yy'])
#bulge:
bulge = dp.getXYdata('../testing/891/891_dtBulge.dat')
b = np.asarray(bulge['yy'])
#gas:
gas = dp.getXYdata('../testing/891/891_dtGas.dat')
g = np.asarray(gas['yy'])

#***************************define total curve
#D=9.25 #disk M-L ratio provided in [1] 
#B=.5 #bulge M-L ratio provided in [1] 

def t(r,B,D,rc,rho00,G):
    return np.sqrt((D*d)**2
                   +(B*b)**2
                   +(nf.h_v(r,rc,rho00))**2
                   +(G*g)**2)


# In[3]:


rc= 1.9#default value
rho00=1e8 #default value

v_err1=v_err1
weighdata=1/v_err1
# LMFit

#Setup
g_mod = lm.Model(t)
g_params = g_mod.make_params()

#gas
g_params.add('G', value=1, min=.99,max=1.01)          #Prefactor
#Bulge
g_params.add('B', value=1, min=0)          #Prefactor

#Disk
g_params.add('D', value=1, min=0)       #Prefactor
#Halo
g_params.add('rc', value=rc, min=0)          #Core radius (kpc)
g_params.add('rho00', value=rho00, min=0)     #Central density 

#Do fit
g_fit = g_mod.fit(v_dat,g_params,r=r,weights=weighdata)


# In[4]:


# Define for plotting
bestg = g_fit.best_fit
#delg = g_fit.eval_uncertainty()
print('Fit information for all-component fit:')
g_fit


# In[5]:
g_dict = g_fit.best_values
g_b = g_dict['B']
g_d = g_dict['D']
g_rc = g_dict['rc']
g_rho00 = g_dict['rho00']
g_g = g_dict['G']
halo_curve = nf.h_v(r,g_rc,g_rho00)