###############
### Imports ###
###############

import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interactive, fixed, FloatSlider, HBox, Layout, Button, Label, Output, VBox
from IPython.display import display
import fitting_NGC891 as fitting891
import NGC5533_functions as nf
import sys
sys.path.append('../../python/')          # Defining path
import dataPython as dp

#############################
### Parameters for slider ###
#############################
      
# Density:
minrho = 0.1                 
maxrho = 3.8                 
defaultrho = 1.5               # default density value for slider
stepRHO = minrho               # step

# For cutoff radius
minrcut = 0.1                 # min number of black holes at the center
maxrcut = 3.0                 # max number of black holes at the center
defaultrcut = 1.4             # default number of bh's at the center for slider
stepRCUT = minrcut            # step of # of bh's at the center for slider

####################################
### Plotting function for widget ###
####################################

def f(bpref,dpref,rcut,rho0,gpref):
    
    # Define r
    r = np.linspace(0.0001,16.4,1000)
    
    # Plot
    plt.figure(figsize=(9,7))
    plt.xlim(0,16)
    plt.ylim(0,255)
    
    plt.errorbar(fitting891.r_dat,fitting891.v_dat,yerr=fitting891.v_err1,fmt='bo',label='Data')
    plt.plot(r,bpref*fitting891.bulge(r),label=("Bulge"),color='orange')
    plt.plot(r,dpref*fitting891.disk(r),label=("Disk"),color='purple')
    plt.plot(r,nf.h_v(r,rcut,rho0),label=("Halo"),color='green')
    plt.plot(r,gpref*fitting891.gas(r),label=("Gas"),color='blue')
    plt.plot(r,fitting891.t(r,bpref,dpref,rcut,rho0,gpref),label=("Total Curve"),color='red')
    plt.title("Interactive Rotation Curve - Galaxy: NGC 891")
    plt.xlabel("Radius (kpc)")
    plt.ylabel("Velocity (km/s)") 
    plt.legend()
    
    # Residuals
    residuals = fitting891.v_dat - fitting891.t(fitting891.r_dat,bpref,dpref,rcut,rho0,gpref)
 
    # Determining errors
    errors = fitting891.v_err1**2 #second term is inclination uncertainty
    # Chi squared
    chisquared = np.sum(residuals**2/errors**2)
    reducedchisquared = chisquared * (1/(len(fitting891.r_dat)-6))
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    plt.text(1,240,r"Reduced $\chi^2$:"+'\n'+"{:.2f}".format(reducedchisquared),ha='left',va='top',bbox=props,size=10)
       
style = {'description_width': 'initial'}
layout = {'width':'800px'}

################################
######## Define Sliders ########
################################

bpref = FloatSlider(min=0, max=5, step=0.05, 
                    value=fitting891.best_bpref, 
                    description='Bulge Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

dpref = FloatSlider(min=0, max=5, step=0.05, 
                    value=fitting891.best_dpref, 
                    description='Disk Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

rcut = fixed(fitting891.best_rcut)

rho0 = FloatSlider(min=1e6, max=5e8, step=1e3, 
                    value=fitting891.best_rho0, 
                    description='Halo Surface Density [$M_{\odot} / pc^3$]', 
                    readout_format='.2e', 
                    orientation='horizontal', 
                    style=style, layout=layout)

gpref = FloatSlider(min=0, max=5, step=0.1, 
                    value=fitting891.best_gpref, 
                    description='Gas Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)


def interactive_plot(f):
    interact = interactive(f,  bpref = bpref, 
                               dpref = dpref, 
                               rcut = rcut,
                               rho0 = rho0,
                               gpref = gpref,
                               continuous_update=False)
    return interact

################################
########### Button #############
################################

# Button to revert back to Best Fit
button = Button(
    description="Best Fit",
    button_style='warning', # 'success', 'info', 'warning', 'danger' or ''
    icon='check')
out = Output()

def on_button_clicked(_):
    bpref.value = fitting891.best_bpref
    dpref.value = fitting891.best_dpref
    rho0.value = fitting891.best_rho0
    gpref.value = fitting891.best_gpref
button.on_click(on_button_clicked)