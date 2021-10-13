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
import ngc891code as wi

# Define components
def bulge(B):
    return B*wi.b

def disk(D):
    return D*wi.d

def halo(r,rc,rho00):
    return nf.h_v(r,rc,rho00,load=False)

def gas(G):
    return G*wi.g

def totalcurve(r,G,B,D,rc,rho00):
    return np.sqrt((gas(G)**2)
                   + (bulge(B)**2) 
                   + (disk(D)**2)
                   + (halo(r,rc,rho00)**2))


# In[9]:


#best fitted prefactor values for each component, to be used as default (initial) values for widget sliders
g_dict = wi.g_fit.best_values
best_G = g_dict['G']
best_B = g_dict['B']
best_D = g_dict['D']
best_rc = g_dict['rc']
best_rho00 = g_dict['rho00']


# In[10]:

z=wi.r[0]
rr=[0,z]

# Define plotting function
def f(G,B,D,rc,rho00):
    
    # Define r
    xmax=20 #kpc
    rval = np.linspace(0,11.2,19)
    
    # Plot
    plt.figure(figsize=(11,6))
    plt.xlim(0,xmax)
    plt.ylim(0,360)
    
    plt.errorbar(wi.r,wi.v_dat,yerr=wi.v_err1,fmt='bo',label='Data')
    
    plt.plot(wi.r,bulge(B),label=("Bulge"),color='orange')
    plt.plot(rr,[0,bulge(B)[0]],color='orange') #straight lining connecting halo curve to origin, for visual aesthetic
    plt.plot(wi.r,disk(D),label=("Disk"),color='purple')
    plt.plot(rr,[0,disk(D)[0]],color='purple') #straight lining connecting halo curve to origin, for visual aesthetic
    plt.plot(wi.r,halo(wi.r,rc,rho00),label=("Halo"),color='green')
    plt.plot(rr,[0,halo(wi.r,rc,rho00)[0]],color='green') #straight lining connecting halo curve to origin, for visual aesthetic
    plt.plot(wi.r,gas(G),label=("Gas"),color='blue')
    plt.plot(wi.r,totalcurve(wi.r,G,B,D,rc,rho00),label=("Total Curve"),color='red')
    plt.plot(rr,[0,totalcurve(wi.r,G,B,D,rc,rho00)[0]],color='red') #straight lining connecting halo curve to origin, for visual aesthetic

    plt.title("Interactive Rotation Curve - Galaxy: NGC7814")
  
    plt.xlabel("Radius (kpc)")
    plt.ylabel("Velocity (km/s)")
    
    
    plt.legend(loc='lower right')
    
     
    #chi squared button commented out for now bc tiny graph
    residuals = wi.v_dat - totalcurve(wi.r,G,B,D,rc,rho00)
    # Determining errors
    errors = wi.v_err1**2 #second term is inclination uncertainty
    # Chi squared
    chisquared = np.sum(residuals**2/errors**2)
    #chisquared = stats.chisquare(v_dat,totalcurve(r,M,bpref,dpref,rc,rho00,gpref))
    reducedchisquared = chisquared * (1/(len(wi.r)-6))
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    plt.text(10,300,r"$\chi^2$: {:.5f}".format(chisquared)+'\n'+r"Reduced: {:.5f}".format(reducedchisquared),bbox=props)
    plt.annotate('Data source: Fraternali1, Sancisi, Kamphuis. "A tale of two galaxies: light and mass in NGC 891 and NGC 7814". A&A Journal. Jun 2011',
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

G = FloatSlider(min=0, max=5, step=0.1, 
                    value=best_G,
                    description='Gas Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

B = FloatSlider(min=0, max=10, step=0.1, 
                    value=best_B, 
                    description='Bulge Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

D = FloatSlider(min=0, max=5, step=0.1, 
                    value=best_D, 
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
                               G = G,
                               B = B, 
                               D = D, 
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
    B.value = best_B
    D.value = best_D
    rho00.value = best_rho00
    G.value = best_G

button.on_click(on_button_clicked)

# displaying button and its output together
VBox([button,out,interactive_plot(f)])


# In[ ]:




