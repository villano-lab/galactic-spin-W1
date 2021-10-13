#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Imports
import sys
sys.path.append('../python/')
#import NGC5533_functions-newmag as nf

import numpy as np
import matplotlib.pyplot as plt
#import scipy.optimize as opt
import lmfit as lm
import dataPython as dp
import scipy.interpolate as inter

import NGC5533_functions as nf

from datetime import datetime
import time
from ipywidgets import interactive, fixed, FloatSlider, HBox, Layout, Button, Label, Output, VBox

from IPython.display import display, clear_output
from IPython.display import Javascript

import scipy.stats as stats
import warnings
warnings.filterwarnings("ignore")  #ignore warnings

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
import widget2code as wi

# Define components
def bulge(BX):
    return BX*wi.rbulgev_fit

def disk(DX):
    return DX*wi.rdiskv_fit

def halo(r_dat,rc,rho00):
    return nf.h_v(r_dat,rc,rho00,load=False)

def gas(GX):
    return GX*wi.rgasv_fit

def totalcurve(rval,GX,BX,DX,rc,rho00):
    return np.sqrt((gas(GX)**2)
                   + (bulge(BX)**2) 
                   + (disk(DX)**2)
                   + (halo(wi.rval,rc,rho00)**2))


# In[9]:


#best fitted prefactor values for each component, to be used as default (initial) values for widget sliders
g_dict = wi.g_fit.best_values
best_GX = g_dict['GX']
best_BX = g_dict['BX']
best_DX = g_dict['DX']
best_rc = g_dict['rc']
best_rho00 = g_dict['rho00']


# In[10]:


# Define plotting function
def f(GX,BX,DX,rc,rho00):
    
    # Define r
    xmax=12 #kpc
    r = np.linspace(0,11.2,19)
    
    # Plot
    plt.figure(figsize=(11,6))
    plt.xlim(0,xmax)
    plt.ylim(0,360)
    
    plt.errorbar(wi.r_dat,wi.v_dat,yerr=wi.v_err1,fmt='bo',label='Data')
    plt.plot(wi.rbulger,bulge(BX),label=("Bulge"),color='orange')
    plt.plot(wi.rdiskr,disk(DX),label=("Disk"),color='purple')
    plt.plot(wi.rval,halo(wi.rval,rc,rho00),label=("Halo"),color='green')
    plt.plot(wi.rgasr,gas(GX),label=("Gas"),color='blue')
    plt.plot(wi.r_dat,totalcurve(wi.rval,GX,BX,DX,rc,rho00),label=("Total Curve"),color='red')
    plt.title("Interactive Rotation Curve - Galaxy: NGC 5005")
    
    #plt.plot(rval,halo_curve,'r-',label='Correct Halo')
    #plt.plot(rbulger,g_b*rbulgev,label='Correct bulge')
    plt.xlabel("Radius (kpc)")
    plt.ylabel("Velocity (km/s)")
    
    
    plt.legend(loc='lower right')
 
    # Chi squared and reduced chi squared
    # Residuals
    
    residuals = wi.v_dat - totalcurve(wi.rval,GX,BX,DX,rc,rho00)
    # Determining errors
    errors = wi.v_err1**2 #second term is inclination uncertainty
    # Chi squared
    chisquared = np.sum(residuals**2/errors**2)
    #chisquared = stats.chisquare(v_dat,totalcurve(r,M,bpref,dpref,rc,rho00,gpref))
    reducedchisquared = chisquared * (1/(len(wi.r_dat)-6))
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    plt.text(10,300,r"$\chi^2$: {:.5f}".format(chisquared)+'\n'+r"Reduced: {:.5f}".format(reducedchisquared),bbox=props)
    plt.annotate('Data Source: Richards, et al. "Baryonic distributions in the dark matter halo of NGC 5005", MNRAS, Volume 449, Issue 4, 01 June 2015, Pages 3981–3996',
            xy=(0, 0), xytext=(0,5),
            xycoords=('axes fraction', 'figure fraction'),
            textcoords='offset points',
            size=10, ha='left', va='bottom')
    
    plt.show()


# In[11]:


# Appearance
style = {'description_width': 'initial'}
layout = {'width':'600px'}

# Define slides

GX = FloatSlider(min=0, max=5, step=0.1, 
                    value=best_GX,
                    description='Gas Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

BX = FloatSlider(min=0, max=5, step=0.1, 
                    value=best_BX, 
                    description='Bulge Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

DX = FloatSlider(min=0, max=5, step=0.1, 
                    value=best_DX, 
                    description='Disk Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

#rc = FloatSlider(min=0, max=5, step=0.1, value=best_rc, description='Halo Core Radius [kpc]', readout_format='.2f', orientation='horizontal', style=style, layout=layout)
rc = fixed(best_rc)
#rho00 = fixed(best_rho00)

rho00 = FloatSlider(min=0, max=1e9, step=1e7, 
                    value=best_rho00, 
                    description=r'Halo Surface Density [$M_{\odot} / pc^3$]', 
                    readout_format='.2e', 
                    orientation='horizontal', 
                    style=style, layout=layout)

# Interactive widget
def interactive_plot(f):
    interact = interactive(f,
                               GX = GX,
                               BX = BX, 
                               DX = DX, 
                               rc = rc,
                               rho00 = rho00,
        
                               continuous_update=False)
    return interact

# Button to revert back to Best Fit
button = Button(
    description="Best Fit",
    button_style='warning', # 'success', 'info', 'warning', 'danger' or ''
    icon='check')
out = Output()

def on_button_clicked(_):
    #display(Javascript('IPython.notebook.execute_cells_below()'))
    BX.value = best_BX
    DX.value = best_DX
    rho00.value = best_rho00
    GX.value = best_GX

button.on_click(on_button_clicked)

# displaying button and its output together
VBox([button,out,interactive_plot(f)])


# In[ ]:




