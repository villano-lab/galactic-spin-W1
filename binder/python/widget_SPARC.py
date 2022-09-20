"""A module for handling the SPARC widget, notebook `09_Widget_SPARC_Galaxies.ipynb <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/09_Widget_SPARC_Galaxies.ipynb>`_.
"""

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

import sys
sys.path.append('python')
import components as comp
px = 1/plt.rcParams['figure.dpi']

################################
############ Data ##############
################################
# Import data from previously downloaded sparc website: http://astroweb.cwru.edu/SPARC/
# Lelli, Federico, Stacy S. McGaugh, and James M. Schombert. “SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves.” The Astronomical Journal 152, no. 6 (2016): 157. https://doi.org/10.3847/0004-6256/152/6/157. 

# Read the file where the chosen galaxy name is located (this was modified by the user)
textfile = open('python/chosengalaxy.txt', 'r')
"""A local text file, opened in read-only mode, indicating a chosen galaxy.

:type: file

The purpose of this file is to:  
1. Allow this chosen galaxy variable to be passed flexibly between programs.
2. Allow this chosen galaxy to be retained between sessions.

.. seealso:: The data in this text file is stored in the :func:`galaxy <widget_SPARC.galaxy>` variable.
"""

galaxy = textfile.read()
"""The name of the chosen galaxy.

:type: string

.. seealso:: This string is retrieved from :func:`textfile <widget_SPARC.textfile>`.
"""

textfile.close()

# Define SPARC directory
#SPARC_file_directory='./data/sparc/'                       #note that '' means the string variable is blank
"""The directory containing SPARC files.

:type: string
"""

# Define file path for .dat files
SPARC_file_path = './data/sparc/' + galaxy + '_rotmod.dat'
"""The filename (including relative path) containing data for the :func:`chosen galaxy <galaxy>`.

:type: string
"""

# Load the galaxy data
data = np.loadtxt(SPARC_file_path)
"""A numpy array of data loaded from the :func:`SPARC file path <SPARC_file_path>`.

:type: array
"""

# Split columns into arrays
Rad,Vobs,errV,Vgas,Vdisk,Vbul,SBdisk,SBbul = data.T

# Check if it has all components
if np.sum(Vbul) == 0:     # Bulge
    warning_bulge = "There is no bulge component."
    """If the bulge component is missing, a message indicating that it is. Otherwise, an empty string.

    :type: string
    """
else: 
    warning_bulge = ""

if np.sum(Vdisk) == 0:    # Disk
    warning_disk = "There is no disk component."
    """If the disk component is missing, a message indicating that it is. Otherwise, an empty string.

    :type: string
    """
else: 
    warning_disk = ""

if np.sum(Vgas) == 0:     # Gas
    warning_gas = "There is no gas component."
    """If the gas component is missing, a message indicating that it is. Otherwise, an empty string.

    :type: string
    """
else: 
    warning_gas = ""

# Define distance to galaxy in Mpc
with open(SPARC_file_path) as file:
    distance = float(file.readline().split()[3])
    """Distance to galaxy, in Mpc.

    :type: float
    """

################################
######### Components ###########
################################

# Bulge
def bulge(r,bpref):
    """
    Interpolating the bulge velocity.

    :parameters:
        r : [array]
            Sampling radius values or distance from the center of the galaxy (in kpc).           
        bpref : [float]
            Bulge prefactor or scaling factor (unitless). 

    :returns:
        Splined bulge velocity as a [function] of sampling radii.

    :example:
        >>> # Define measured radius and velocity and interpolate them
        >>> import numpy as np
        >>> Rad = np.array([0, 0.5, 1.2, 2.6, 5.3, 6.7, 7.1, 9.5])
        >>> Vbul = np.array([100, 200, 280, 200, 130, 120, 110, 100])
        >>> r = np.linspace(0, 20, 50)
        >>> bulgespline = bulge(r, bpref=1)
        >>> plt.plot(Rad, Vbul, 'ro')
        >>> plt.plot(r, bulgespline)
        >>> plt.show()

        .. image:: ../../images/bulgeExample.png
            :alt: A plot of the supplied bulge velocities, with radii as the x-points and velocities as the y-points, and a blue line representing a rough interpolation of them with scaling factor 1.
    """
    
    polynomial = InterpolatedUnivariateSpline(Rad,bpref*Vbul,k=3)
    return polynomial(r)

# Disk
def disk(r,dpref):
    """
    Interpolating the disk velocity.

    :parameters:
        r : [array]
            Sampling radius values or distance from the center of the galaxy (in kpc).           
        dpref : [float]
            Disk prefactor or scaling factor (unitless). 

    :returns:
        Splined disk velocity as a [function] of sampling radii.

    :example:
        >>> # Define measured radius and velocity and interpolate them
        >>> #Rad = #Plotting Rad and Vdisk as imported from this library
        >>> #Vdisk = #^
        >>> r = np.linspace(0, 20, 50)
        >>> diskspline = disk(r, dpref=1)
        >>> plt.plot(Rad, Vdisk, 'ro')
        >>> plt.plot(r, diskspline)
        >>> plt.show()

        .. image:: ../../images/diskExample.png
            :alt: A plot of the library's default disk velocities, with radii as the x-points and velocities as the y-points, and a blue line representing a somewhat rough interpolation of them with scaling factor 1.
    """
    
    polynomial = InterpolatedUnivariateSpline(Rad,dpref*Vdisk,k=3)   
    return polynomial(r)

# Gas 
def gas(r):
    """
    Interpolating the gas velocity.

    :parameters:
        r : [array]
            Sampling radius values or distance from the center of the galaxy (in kpc).

    :returns:
        Splined gas velocity as a [function] of sampling radii.

    :example:
        >>> # Define measured radius and velocity and interpolate them
        >>> #Rad = # We plot the radius and velocities as imported from this library.
        >>> #Vgas = #^
        >>> r = np.linspace(0, 20, 50)
        >>> gasspline = gas(r)
        >>> plt.plot(Rad, Vgas, 'ro')
        >>> plt.plot(r, gasspline)
        >>> plt.show()

        .. image:: ../../images/gasExample.png
            :alt: A plot of the library's default disk velocities, with radii as the x-points and velocities as the y-points, and a blue line representing a somewhat rough interpolation of them with scaling factor 1.
    """
    
    polynomial = InterpolatedUnivariateSpline(Rad,Vgas,k=3)   
    return polynomial(r)

# Halo

## Equation for Dark Matter halo velocity
def halo(r,
         rc,
         rho0):
    """
    Function to calculate the gravitational effect of a Dark Matter halo using the isothermal density profile (Source: [Jimenez2003]_).

    :parameters:
        r : [array]
            Radius values or distance from the center of the galaxy used to calculate velocities (in kpc). 
        rc : [float]
            Cutoff radius (in kpc). Default: `1.4`
        rho0 : [float]
            Central mass density (in solar mass/kpc^3). Default: `0.31e9`

    :returns:
        A [float] or an [array] of halo velocities (in km/s).

    :example:
        >>> # Calculate the gravitational effect of the Dark Matter halo 
        >>> # of NGC 5533, 10 kpc away. 
        >>> print(halo(r=np.array([10,15,20,25,30,35,40,45,50,100]), 
        ...     rc=1.4, rho0=0.31e9)))
        [162.0220417  168.23695403 171.41313542 173.33987823 174.63289949 
        175.5605844  176.25855891 176.80272454 177.23886723 179.21029129]
    """ 
    
    return np.sqrt(4*np.pi*comp.G*rho0*(rc**2)*(1-((rc/r)*np.arctan(r/rc))))

# Total velocity
def totalcurve(r,
               bpref,
               dpref,
               rc,
               rho0):
    """
    Function to calculate the total gravitational effect of all components of a galaxy. 
    The velocities of each component is added in quadrature to calculate the total rotational velocity.

    :parameters:
        r : [array]
            Radius values or distance from the center of the galaxy used to calculate velocities (in kpc).            
        bpref : [float]
            Bulge prefactor or scaling factor (unitless).       
        dpref : [float]
            Disk prefactor or scaling factor (unitless). 
        rc : [float]
            Cutoff radius (in kpc).          
        rho0 : [float]
            Central mass density (in solar mass/kpc^3).

    :returns:
        A [float] or an [array] of total velocities (in km/s).

    .. note::
        If the galaxy contains a supermassive black hole at the center, it is incorporated in the bulge velocity.

    :example:
        >>> # Calculate the gravitational effect of all components of a galaxy 
        >>> #at the distance of 10,15,20,25,30,35,40,45,50, and 100 kpc. 
        >>> print(totalcurve(r=np.array([10,15,20,25,30,35,40,45,50,100]), 
        ...     bpref=1, dpref=1, rc=1.4, rho0=0.31e9)))
        [3.34371479e+02 3.22215072e+02 7.25902496e+02 2.16917607e+03 5.21519130e+03 
        1.04268198e+04 1.83732304e+04 2.96250118e+04 4.47531137e+04 5.34875631e+05]
    """
    
    # Total velocity with components added in quadrature
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
"""An lmfit Model object for the :func:`totalcurve <widget_SPARC.totalcurve>` function.

:type: lmfit.Model
"""
fit_params = fit_mod.make_params()
"""An lmfit Parameters object for the :func:`fit_mod <widget_SPARC.fit_mod>` Model.

:type: lmfit.Parameters

+----------+----------------+---------+---------+
| Parameter| Starting Value | Minimum | Maximum |
+==========+================+=========+=========+
| rc       | 1.4            | 0.1     | (None)  |
+----------+----------------+---------+---------+
| rho0     | 3.10e8         | 0       | (None)  |
+----------+----------------+---------+---------+
| bpref    | 1              | 0.5     | 100     |
+----------+----------------+---------+---------+
| dpref    | 1              | 0.5     | 100     |
+----------+----------------+---------+---------+

"""
#weighdata = 1/errV

# Halo
fit_params.add('rc',    value=1.4,    min=0.1)           # Cutoff radius (kpc)
fit_params.add('rho0',  value=3.10e8,  min=0)             # Central halo mass density (Solar mass/kpc^3)
# Bulge
fit_params.add('bpref', value=1,     min=0.5, max=100)  # Bulge prefactor
# Disk
fit_params.add('dpref', value=1,     min=0.5, max=100)  # Disk prefactor

# Do fit
fit = fit_mod.fit(Vobs,fit_params,r=Rad,weights=1/errV)
"""An lmfit fitting result for the :func:`fit_mod <widget_SPARC.fit_mod>` Model against the :func:`data imported for the chosen galaxy <widget_SPARC.data>`.

:type: lmfit.ModelResult
"""

#################################
### Define fitting parameters ###
#################################

bestfit = fit.best_fit
"""Best fit results of :func:`the fit of the total curve against the galaxy's data <widget_SPARC.fit>`.

:type: ndarray
"""

fit_dict = fit.best_values
"""Dictionary of best parameter results from :func:`the fit of the total curve against the galaxy's data <widget_SPARC.fit>`.

:type: dict 

:keys: `bpref`, `dpref`, `rc`, `rho0`.
"""

#####################
###### Widget #######
#####################

# Define plotting function
def widgetfunction(bpref,dpref,rc,rho0):
    """
    Generate a plot for use with interactive rotation curve plot of SPARC data with sliders for parameters.
    Can also be used to generate a static plot of individual components and their total.
    
    :parameters:
        bpref: [float]
            Prefactor scaling the bulge component.
        dpref: [float]
            Prefactor scaling the disk component.
        rcut: [float]
            Cutoff radius of the halo.
        rho0: [float]
            Density parameter for the halo.

    :returns:
        None

    .. seealso:: For information on how the curves displayed are calculated, see: 
        :func:`totalcurve <widget_SPARC.totalcurve>`, 
        :func:`bulge <widget_SPARC.bulge>`,
        :func:`disk <widget_SPARC.disk>`,
        :func:`halo <widget_SPARC.halo>`,
        :func:`gas <widget_SPARC.gas>`.  

        See the :func:`interactive_plot <widget_SPARC.interactive_plot>` function as an example usecase of this function.
    """
    
    # Define radius
    r = np.linspace(np.min(Rad),np.max(Rad),1000)
    
    # Plot
    plt.figure(figsize=(1100*px,700*px))
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
            xy=(0, 0), xytext=(0,0),
            xycoords=('axes fraction', 'figure fraction'),
            textcoords='offset points',
            size=12, ha='left', va='bottom')
    plt.tight_layout()
    plt.legend(bbox_to_anchor=(1,1), loc="upper left") 
    plt.show()

################################
######### Interactive ##########
################################

bpref = FloatSlider(min=0, max=5, step=0.1, 
                    value=fit_dict['bpref'], 
                    description='Bulge Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style={'description_width': 'initial'}, 
                    layout={'width':'600px'},disabled=not bool(np.sum(Vbul)))
"""Slider for controlling the bulge prefactor.

:type: ipywidgets.widgets.widget_float.FloatSlider
"""

dpref = FloatSlider(min=0, max=5, step=0.1, 
                value=fit_dict['dpref'], 
                description='Disk Prefactor', 
                readout_format='.2f', 
                orientation='horizontal', 
                style={'description_width': 'initial'}, 
                layout={'width':'600px'},disabled=not bool(np.sum(Vdisk)))
"""Slider for controlling the disk prefactor.

:type: ipywidgets.widgets.widget_float.FloatSlider
"""

rc = FloatSlider(min=0, max=20, step=0.1, 
                value=fit_dict['rc'], 
                description='Halo Core Radius [$kpc$]', 
                readout_format='.2f', 
                orientation='horizontal', 
                style={'description_width': 'initial'}, layout={'width':'600px'})
"""Slider for controlling the cutoff radius.

:type: ipywidgets.widgets.widget_float.FloatSlider
"""

rho0 = FloatSlider(min=0, max=1e11, step=1e5, 
                value=fit_dict['rho0'], 
                description='Halo Central Mass Density [$M_{\odot} / kpc^3$]', 
                readout_format='.2e', 
                orientation='horizontal', 
                style={'description_width': 'initial'}, layout={'width':'600px'})
"""Slider for controlling the halo density variable rho0.

:type: ipywidgets.widgets.widget_float.FloatSlider
"""

def interactive_plot(widgetfunction):
    """
    Generate an interactive plot widget, allowing the user to interact with the SPARC data and the galaxy's components.

    :parameters:
        widgetfunction: [function]
            A function that generates the base plot for the widget to alter. This should, in all likelihood, be :func:`widgetfunction <widget_SPARC.widgetfunction>`.

    :returns: 
        [ipywidgets.widgets.interaction.interactive] -- creates sliders to make the plot interactive.

    .. seealso:: For an example usage of this function, see the notebook `09_Widget_SPARC_Galaxies.ipynb on Binder <https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/HEAD?labpath=binder%2F09_Widget_SPARC_Galaxies.ipynb>`_.

    """

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
"""A button that returns all settings to the best fit.

:type: ipywidgets.widgets.widget_button.Button
"""
out = Output()
"""A handler for widget output.

:type: ipywidgets.widgets.widget_output.Output
"""

def on_button_clicked(_):
    """
    A function to reset values when the 'Best Fit' button is clicked. 
    
    :parameters: None

    :returns:
        None

    :example:
        >>> button.on_click(on_button_clicked)

        This renders the button click behavior seen in `09_Widget_SPARC_Galaxies.ipynb on Binder <https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/HEAD?labpath=binder%2F09_Widget_SPARC_Galaxies.ipynb>`__.
    """
    bpref.value = fit_dict['bpref']
    dpref.value = fit_dict['dpref']
    rc.value = fit_dict['rc']
    rho0.value = fit_dict['rho0']

button.on_click(on_button_clicked)


################################
######## Galaxy Image ##########
################################

from astroquery.skyview import SkyView
from astropy.wcs import WCS

def GalaxyImage(position=galaxy,survey=['DSS']):
    """
    Fetch and display an image of the selected galaxy.

    :parameters:
        Position: [string]
            The name of the galaxy whose image is to be fetched.
        Survey: [list]
            The name(s) of surveys, as strings, containing the galaxy whose image is to be fetched.

    :returns:
        None

    .. seealso:: For an example usage of this function, see the notebook `09_Widget_SPARC_Galaxies.ipynb on Binder <https://mybinder.org/v2/gh/villano-lab/galactic-spin-W1/HEAD?labpath=binder%2F09_Widget_SPARC_Galaxies.ipynb>`__.
    """
    # DSS images of the target
    hdu = SkyView.get_images(galaxy, survey=survey)[0][0]
    gfilter = hdu.data

    # WCS 
    plt.figure(figsize=(500*px,500*px))
    wcs = WCS(hdu.header)
    ax = plt.gca(projection=wcs)

    # Plot galaxy image
    ax.imshow(hdu.data, vmin=np.percentile(gfilter,0), vmax=np.percentile(gfilter,100), cmap='plasma')
    ax.set(xlabel="RA", ylabel="Dec")
    plt.title("{}".format(galaxy),fontsize='14')
    plt.show()