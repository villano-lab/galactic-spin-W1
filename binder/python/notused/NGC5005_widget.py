###############
### Imports ###
###############
import numpy as np
import matplotlib.pyplot as plt
import lmfit as lm
import dataPython as dp
import scipy.interpolate as inter
import NGC5005_components as f

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

def totalcurve(r, bpref,dpref,gpref,rc,rho0):
    return np.sqrt((f.gas(r,gpref)**2)
                   + (f.bulge(r,bpref)**2) 
                   + (f.disk(r,dpref)**2)
                   + (f.halo(r,rc,rho0)**2))


# Define plotting function
def f(bpref,dpref,gpref,rc,rho00):
    
    # Define r
    xmax=12 #kpc
    r = np.linspace(0,11.2,19)
    
    # Plot
    plt.figure(figsize=(11,6))
    plt.xlim(0,xmax)
    plt.ylim(0,360)
    
    plt.errorbar(fitting.r_dat,fitting.v_dat,yerr=fitting.v_err1,fmt='bo',label='Data')
    plt.plot(r,bulge(r,bpref),label=("Bulge"),color='orange')
    plt.plot(r,disk(r,dpref),label=("Disk"),color='purple')
    plt.plot(r,halo(r,rc,rho00),label=("Halo"),color='green')
    plt.plot(r,gas(r,gpref),label=("Gas"),color='blue')
    plt.plot(r,totalcurve(r,bpref,dpref,gpref,rc,rho0),label=("Total Curve"),color='red')
    plt.title("Interactive Rotation Curve - Galaxy: NGC 5005")
    
    plt.xlabel("Radius (kpc)")
    plt.ylabel("Velocity (km/s)")
    
    plt.legend(loc='lower right')
 
    # Chi squared and reduced chi squared
    # Residuals
    residuals = fitting.v_dat - totalcurve(fitting.rval,bpref,dpref,gpref,rc,rho0)
    # Determining errors
    errors = fitting.v_err1**2 #second term is inclination uncertainty
    # Chi squared
    chisquared = np.sum(residuals**2/errors**2)
    #chisquared = stats.chisquare(v_dat,totalcurve(r,M,bpref,dpref,rc,rho00,gpref))
    reducedchisquared = chisquared * (1/(len(fitting.r_dat)-6))
    
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
gpref = FloatSlider(min=0, max=5, step=0.1, 
                    value=f.gpref,
                    description='Gas Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

bpref = FloatSlider(min=0, max=5, step=0.1, 
                    value=f.bpref, 
                    description='Bulge Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

dpref = FloatSlider(min=0, max=5, step=0.1, 
                    value=f.dpref, 
                    description='Disk Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

rc = FloatSlider(min=0, max=5, step=0.1, 
                 value=f.rc, 
                 description='Halo Core Radius [kpc]', 
                 readout_format='.2f', 
                 orientation='horizontal', 
                 style=style, layout=layout)
#rc = fixed(fitting.rc)
#rho00 = fixed(best_rho00)

rho00 = FloatSlider(min=0, max=1e9, step=1e7, 
                    value=f.rho00, 
                    description=r'Halo Surface Density [$M_{\odot} / pc^3$]', 
                    readout_format='.2e', 
                    orientation='horizontal', 
                    style=style, layout=layout)

# Interactive widget
def interactive_plot(f):
    interact = interactive(f,
                               bpref=bpref,
                               dpref=dpref, 
                               gpref=gpref, 
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
    bpref.value = f.bpref
    dpref.value = f.dpref
    gpref.value = f.gpref
    rc.value=f.rc
    rho00.value = f.rho00

button.on_click(on_button_clicked)