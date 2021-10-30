
#Imports
import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm
import dataPython as dp
import scipy.interpolate as inter
import NGC5533_functions as nf
import NGC5005_components as comp5005

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

# Define components
def bulge(bpref):
    return BX*comp5005.rbulgev_fit

def disk(dpref):
    return DX*comp5005.rdiskv_fit

def halo(r_dat,rc,rho00):
    return nf.h_v(r_dat,rc,rho0,load=False)

def gas(gpref):
    return GX*comp5005.rgasv_fit

def totalcurve(r,bpref,dpref,gpref,rc,rho0):
    return np.sqrt((gas(gpref)**2)
                   + (bulge(bpref)**2) 
                   + (disk(dpref)**2)
                   + (halo(r,rc,rho0)**2))

# Prefactors
bpref = 0.28228633
dpref = 0.71035496
gpref = 0.956

# Define plotting function
def f(gpref,bpref,dpref,rc,rho00):
    
    # Define r
    xmax=12 #kpc
    r = np.linspace(0,11.2,19)
    
    # Plot
    plt.figure(figsize=(11,6))
    plt.xlim(0,xmax)
    plt.ylim(0,360)
    
    plt.errorbar(comp5005.r_dat,comp5005.v_dat,yerr=wi.v_err1,fmt='bo',label='Data')
    plt.plot(comp5005.rbulger,bulge(bpref),label=("Bulge"),color='orange')
    plt.plot(comp5005.rdiskr,disk(dpref),label=("Disk"),color='purple')
    plt.plot(comp5005.rval,halo(comp5005.rval,rc,rho00),label=("Halo"),color='green')
    plt.plot(comp5005.rgasr,gas(gpref),label=("Gas"),color='blue')
    plt.plot(comp5005.r_dat,totalcurve(comp5005.rval,bpref,dpref,gpref,rc,rho0),label=("Total Curve"),color='red')
    plt.title("Interactive Rotation Curve - Galaxy: NGC 5005")
    
    plt.xlabel("Radius (kpc)")
    plt.ylabel("Velocity (km/s)")
    
    plt.legend(loc='lower right')
 
    # Chi squared and reduced chi squared
    # Residuals
    residuals = comp5005.v_dat - totalcurve(comp5005.rval,bpref,dpref,gpref,rc,rho0)
    # Determining errors
    errors = comp5005.v_err1**2 #second term is inclination uncertainty
    # Chi squared
    chisquared = np.sum(residuals**2/errors**2)
    #chisquared = stats.chisquare(v_dat,totalcurve(r,M,bpref,dpref,rc,rho00,gpref))
    reducedchisquared = chisquared * (1/(len(comp5005.r_dat)-6))
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    plt.text(10,300,r"$\chi^2$: {:.5f}".format(chisquared)+'\n'+r"Reduced: {:.5f}".format(reducedchisquared),bbox=props)
    plt.annotate('Data Source: Richards, et al. "Baryonic distributions in the dark matter halo of NGC 5005", MNRAS, Volume 449, Issue 4, 01 June 2015, Pages 3981–3996',
            xy=(0, 0), xytext=(0,5),
            xycoords=('axes fraction', 'figure fraction'),
            textcoords='offset points',
            size=10, ha='left', va='bottom')
    
    plt.show()

# Appearance
style = {'description_width': 'initial'}
layout = {'width':'600px'}

# Define slides
GX = FloatSlider(min=0, max=5, step=0.1, 
                    value=gpref,
                    description='Gas Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

BX = FloatSlider(min=0, max=5, step=0.1, 
                    value=bpref, 
                    description='Bulge Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

DX = FloatSlider(min=0, max=5, step=0.1, 
                    value=dpref, 
                    description='Disk Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

#rc = FloatSlider(min=0, max=5, step=0.1, value=best_rc, description='Halo Core Radius [kpc]', readout_format='.2f', orientation='horizontal', style=style, layout=layout)
rc = fixed(comp5005.rc)
#rho00 = fixed(best_rho00)

rho0 = FloatSlider(min=0, max=1e9, step=1e7, 
                    value=comp5005.rho0, 
                    description=r'Halo Surface Density [$M_{\odot} / pc^3$]', 
                    readout_format='.2e', 
                    orientation='horizontal', 
                    style=style, layout=layout)

# Interactive widget
def interactive_plot(f):
    interact = interactive(f,
                               bpref,
                               dpref, 
                               gpref, 
                               rc = rc,
                               rho0 = rho0,
        
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
    BX.value = bpref
    DX.value = dpref
    rho00.value = best_rho00
    GX.value = best_GX

button.on_click(on_button_clicked)

# displaying button and its output together
VBox([button,out,interactive_plot(f)])