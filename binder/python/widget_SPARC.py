################################
########### Imports ############
################################

import numpy as np                                              # Basic operators
import matplotlib.pyplot as plt                                 # Plotting
import os                                                       # Listing filenames

from scipy.interpolate import InterpolatedUnivariateSpline      # Interpolation (spline) function
import lmfit as lm                                              # Fitting

from ipywidgets import interactive, fixed, FloatSlider, HBox, Layout, Button, Label, Output, VBox
from IPython.display import display

import warnings
warnings.filterwarnings("ignore")         # ignore warnings

################################
############ Data ##############
################################
# Import data from previously downloaded sparc website: http://astroweb.cwru.edu/SPARC/
# Lelli, Federico, Stacy S. McGaugh, and James M. Schombert. “SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves.” The Astronomical Journal 152, no. 6 (2016): 157. https://doi.org/10.3847/0004-6256/152/6/157. 

# Read the file where the chosen galaxy name is located (this was modified by the user)
textfile = open('python/chosengalaxy.txt', 'r')
galaxy = textfile.read()
textfile.close()

# Define SPARC directory
SPARC_file_directory='./data/sparc/'                       #note that '' means the string variable is blank

# Define file path for .dat files
SPARC_file_path = SPARC_file_directory + galaxy + '_rotmod.dat'

# Load the galaxy data
data = np.loadtxt(SPARC_file_path)

# Split columns into arrays
Rad,Vobs,errV,Vgas,Vdisk,Vbul,SBdisk,SBbul = data.T

# Check if it has all components
if np.sum(Vbul) == 0:     # Bulge
    nobulge = True
    warning_bulge = "There is no bulge component."
else: 
    nobulge = False
    warning_bulge = ""

if np.sum(Vdisk) == 0:    # Disk
    nodisk = True
    warning_disk = "There is no disk component."
else: 
    nodisk = False
    warning_disk = ""

if np.sum(Vgas) == 0:     # Gas
    nogas = True
    warning_gas = "There is no gas component."
else: 
    nogas = False
    warning_gas = ""

# Define distance to galaxy in Mpc
firstline = open(SPARC_file_path).readline()
firstline = firstline.split()
distance = float(firstline[3])

##############################
### Interpolation function ###
##############################

def interpd(x,y):
    return InterpolatedUnivariateSpline(x,y,k=3)

################################
######### Components ###########
################################

# Bulge
def bulge(r,bpref):
    polynomial = interpd(Rad,bpref*Vbul)   
    return polynomial(r)

# Disk
def disk(r,dpref):
    polynomial = interpd(Rad,dpref*Vdisk)   
    return polynomial(r)

# Gas 
def gas(r):
    polynomial = interpd(Rad,Vgas)   
    return polynomial(r)

# Halo 
## Set parameters
rho0 = 3.10e8       # Central mass density (in solar mass/kpc^3)
rc = 1.4            # Core radius (in kpc)
G = 4.300e-6        # Gravitational constant (kpc/solar mass*(km/s)^2)

## Equation for dark matter halo velocity
def halo(r,rc,rho0):
    v = np.sqrt(4*np.pi*G*rho0*rc**2*(1 - rc/r * np.arctan(r/rc)))
    return v

# Total velocity
def totalcurve(r,bpref,dpref,rc,rho0):
    total = np.sqrt(bulge(r,bpref)**2 
                    + disk(r,dpref)**2
                    + halo(r,rc,rho0)**2
                    + gas(r)**2)
    return total

#################################
#### Find Fitting Parameters ####
#################################

# Setup
fit_mod = lm.Model(totalcurve)
fit_params = fit_mod.make_params()
weighdata = 1/errV

# Halo
fit_params.add('rc',    value=rc,    min=0.1)           # Core Radius (kpc)
fit_params.add('rho0',  value=rho0,  min=0)             # Halo Density (Solar mass/kpc^3)
# Bulge
fit_params.add('bpref', value=1,     min=0.5, max=100)  # Bulge Prefactor
# Disk
fit_params.add('dpref', value=1,     min=0.5, max=100)  # Disk Prefactor

# Do fit
fit = fit_mod.fit(Vobs,fit_params,r=Rad,weights=weighdata)

#################################
### Define fitting parameters ###
#################################

bestfit = fit.best_fit

fit_dict = fit.best_values
best_bpref = fit_dict['bpref']
best_dpref = fit_dict['dpref']
best_rc = fit_dict['rc']
best_rho0 = fit_dict['rho0']

#####################
###### Widget #######
#####################

# Define plotting function
def widgetfunction(bpref,dpref,rc,rho0):
    
    # Define radius
    r = np.linspace(np.min(Rad),np.max(Rad),1000)
    
    # Plot
    plt.figure(figsize=(11,7))
    plt.xlim(0,np.max(Rad)+0.2)
    plt.ylim(0,np.max(Vobs)+100)
    
    plt.errorbar(Rad,Vobs,yerr=errV,fmt='bo',label='Data')
    plt.plot(r,bulge(r,bpref),label=("Bulge"),color='orange')
    plt.plot(r,disk(r,dpref),label=("Disk"),color='purple')
    plt.plot(r,halo(r,rc,rho0),label=("Dark Matter Halo"),color='green')
    plt.plot(r,gas(r),label=("Gas"),color='blue')
    plt.plot(r,totalcurve(r,bpref,dpref,rc,rho0),label=("Total Curve"),color='red')
    plt.suptitle("Interactive Rotation Curve - Galaxy: {}".format(galaxy),fontsize='16')
    plt.title("Distance: {} Mpc".format(distance),fontsize='13')
    plt.xlabel("Radius (kpc)",fontsize='14')
    plt.ylabel("Velocity (km/s)",fontsize='14')
    
    # Chi squared and reduced chi squared
    # Residuals
    residuals = Vobs - totalcurve(Rad,bpref,dpref,rc,rho0)
    # Chi squared
    chisquared = np.sum(residuals**2/errV**2)
    dof = len(Rad) - 4                 # number of degrees of freedom = number of observed data - number of fitting parameters
    reducedchisquared = chisquared / dof
    
    if len(Rad) <= 4:
        warning_chi = "Number of measured points are too low to calculate the reduced chi-squared value."
    else:
        warning_chi = ""
  
    plt.annotate(r"$\chi^2$: {:.5f}".format(chisquared)+'\n'
                 +r"Reduced $\chi^2$: {:.5f}".format(reducedchisquared)+'\n'
                 +r"Note: {}".format(warning_chi) 
                 + "{}".format(warning_bulge) 
                 + "{}".format(warning_disk) 
                 + "{}".format(warning_gas),
            xy=(0, 0), xytext=(0,-2),
            xycoords=('axes fraction', 'figure fraction'),
            textcoords='offset points',
            size=13, ha='left', va='bottom')
    plt.legend(bbox_to_anchor=(1,1), loc="upper left") 
    plt.show()

###############################
######### Appearance ##########
###############################

style = {'description_width': 'initial'}
layout = {'width':'600px'}

################################
######## Define Sliders ########
################################

if nobulge == False:  
    bpref = FloatSlider(min=0, max=5, step=0.1, 
                        value=best_bpref, 
                        description='Bulge Prefactor', 
                        readout_format='.2f', 
                        orientation='horizontal', 
                        style=style, layout=layout)
if nobulge == True:
    bpref = FloatSlider(min=0, max=5, step=0.1, 
                        value=best_bpref, 
                        description='Bulge Prefactor', 
                        readout_format='.2f', 
                        orientation='horizontal', 
                        style=style, layout=layout, 
                        disabled=True)

if nodisk == False:
    dpref = FloatSlider(min=0, max=5, step=0.1, 
                        value=best_dpref, 
                        description='Disk Prefactor', 
                        readout_format='.2f', 
                        orientation='horizontal', 
                        style=style, layout=layout)

if nodisk == True:
    dpref = FloatSlider(min=0, max=5, step=0.1, 
                        value=best_dpref, 
                        description='Disk Prefactor', 
                        readout_format='.2f', 
                        orientation='horizontal', 
                        style=style, layout=layout, 
                        disabled=True)    

rc = FloatSlider(min=0, max=20, step=0.1, 
                 value=best_rc, 
                 description='Halo Core Radius [$kpc$]', 
                 readout_format='.2f', 
                 orientation='horizontal', 
                 style=style, layout=layout)
#rc = fixed(best_rc)

rho0 = FloatSlider(min=0, max=1e11, step=1e6, 
                    value=best_rho0, 
                    description='Halo Central Mass Density [$M_{\odot} / kpc^3$]', 
                    readout_format='.2e', 
                    orientation='horizontal', 
                    style=style, layout=layout)


################################
######### Interactive ##########
################################

def interactive_plot(widgetfunction):
    interact = interactive(widgetfunction,  bpref = bpref, 
                               dpref = dpref, 
                               rc = rc,
                               rho0 = rho0,
                               continuous_update=False)
    return interact

################################
########### Button #############
################################

button = Button(
    description="Best Fit",
    button_style='warning', # 'success', 'info', 'warning', 'danger' or ''
    icon='check')
out = Output()

def on_button_clicked(_):
    bpref.value = best_bpref
    dpref.value = best_dpref
    rc.value = best_rc
    rho0.value = best_rho0

button.on_click(on_button_clicked)


################################
######## Galaxy Image ##########
################################

from astroquery.skyview import SkyView
from astropy.wcs import WCS

def GalaxyImage(galaxy,survey='DSS'):
    
    # DSS images of the target
    hdu = SkyView.get_images(galaxy, survey=survey)[0][0]
    gfilter = hdu.data

    # WCS 
    plt.figure(figsize=(5,5))
    wcs = WCS(hdu.header)
    ax = plt.gca(projection=wcs)

    # Plot galaxy image
    ax.imshow(hdu.data, vmin=np.percentile(gfilter,0), vmax=np.percentile(gfilter,100), cmap='plasma')
    ax.set(xlabel="RA", ylabel="Dec")
    plt.title("{}".format(galaxy),fontsize='14')
    plt.show()