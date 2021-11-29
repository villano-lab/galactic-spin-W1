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
import scipy.integrate as si
from scipy.interpolate import InterpolatedUnivariateSpline
import NGC5533_functions as nf             # Components


import scipy.stats as stats
import warnings
warnings.filterwarnings("ignore")  #ignore warnings



# this cell is a placeholder for the .py file

# Read the file where the chosen galaxy name is located
f = open('python/wgalaxy.txt', 'r')
galaxy = f.read()

##############################
### Import data/text files ###
##############################
if galaxy == 'NGC5005':
    dataFile = 'ngc5005_data.txt'
    diskFile='ngc5005_disk.txt'
    bulgeFile = 'ngc5005_bulge.txt'
    gasFile='ngc5005_gas.txt'
    
    dprefs = 0.71035496
    bprefs = 0.28228633
    gprefs = 0.956
    rho0s = 1.5981e+08 #ALSO FROM FITTING?
    rcs = 2.5                  # cutoff radius (in kpc) (source: in paragraph right above fig. 9)

    
    DataSource= 'Richards, et al. "Baryonic distributions in the dark matter halo of NGC 5005", MNRAS, Volume 449, Issue 4, 01 June 2015, Pages 3981–3996'
elif galaxy=='NGC7814':
    dataFile = '7814_data.txt'
    diskFile='7814_disk_gipsy.txt'
    bulgeFile = '7814_bulge_gipsy.txt'
    gasFile='7814_gas_gipsy.txt'
    
    dprefs = 3.85244293       # disk
    bprefs = 4.98331928       # bulge
    gprefs = 1.00239573       # gas
    
    # NGC 7814 (Source: https://www.aanda.org/articles/aa/abs/2011/07/aa16634-11/aa16634-11.html)
    rcs = 2.1                # cutoff radius (in kpc) from Table 5.
    rho0s = 152.3e-3           # central density (in solar mass/pc^3) from Table 5.
    
    DataSource= 'Fraternali1, Sancisi, Kamphuis. "A tale of two galaxies: light and mass in NGC 891 and NGC 7814". A&A Journal. Jun 2011'
elif galaxy=='NGC891':
    dataFile = '891_data.txt'
    diskFile='891_disk_gipsy.txt'
    bulgeFile = '891_bulge_gipsy.txt'
    gasFile='891_gas_gipsy.txt'
    
    #from fitting:
    gprefs=1
    dprefs=0.89567284
    bprefs=0.99506827 
    rho0s = 47435774.7           # this is the value attained from fitting. below is the value provided in paper
    # NGC 891 (Source: https://www.aanda.org/articles/aa/abs/2011/07/aa16634-11/aa16634-11.html)
    rcs=1.9 
    #rho0 = 1e9**33.1e-3        # central density (converted to solar mass/kpc^3) from Table 5.
  
    DataSource= 'Fraternali1, Sancisi, Kamphuis. "A tale of two galaxies: light and mass in NGC 891 and NGC 7814". A&A Journal. Jun 2011'
elif galaxy=='NGC5533':
    dataFile = '100kpc_data.txt' #this or the .txt called 'noord-120kpc-total.txt'???
    diskFile='noord-120kpc-disk.txt'
    bulgeFile = 'noord-120kpc-bulge.txt'
    gasFile='noord-120kpc-gas.txt'
    
    #from fitting:
    gprefs=1
    dprefs=.99
    bprefs=1.03
    rho0s =2.94e9
    #from noord paper:
    #rho0 = .31e9           # central density (converted to solar mass/kpc^3) from Table 5.
    rcs=1.4  
    
    DataSource='E. Noordermeer. The rotation curves of flattened Sérsic bulges. MNRAS,385(3):1359–1364, Apr 2008'
#*******************************
#Assigning txt file info
#*******************************
dataPath=str('data/'+galaxy+'/'+dataFile)
diskPath=str('data/'+galaxy+'/'+diskFile)
bulgePath=str('data/'+galaxy+'/'+bulgeFile)
gasPath=str('data/'+galaxy+'/'+gasFile)

# Datapoints:
data = dp.getXYdata_wXYerr(dataPath)
r_dat = np.asarray(data['xx'])
v_dat = np.asarray(data['yy'])
v_err1 = np.asarray(data['ey'])
# Disk:
disk_data = dp.getXYdata(diskPath)
disk_raw = np.asarray(disk_data['yy'])
# Bulge:
bulge_data = dp.getXYdata(bulgePath)
bulge_raw = np.asarray(bulge_data['yy'])
# Gas:
gas_data = dp.getXYdata(gasPath)
gas_raw = np.asarray(gas_data['yy'])  


print(len(r_dat))
print(len(gas_raw))
print(len(bulge_raw))
print(len(disk_raw))
weighdata = 1/v_err1
G = 4.30091e-6            # gravitational constant (kpc/solar mass*(km/s)^2)

def interpd(x,y):
    return InterpolatedUnivariateSpline(x,y,k=5)

def bulge(r,bpref):
    polynomial = interpd(r_dat,bpref*bulge_raw)   
    return polynomial(r)

def disk(r,dpref):
    polynomial = interpd(r_dat,dpref*disk_raw)   
    return polynomial(r)

def gas(r,gpref):
    polynomial = interpd(r_dat,gpref*gas_raw)   
    return polynomial(r)

def halo(r,rc,rho0):
    return nf.h_v(r,rc,rho0,load=False)


##################################
### Calculating total velocity ###
##################################

def totalvelocity(r,bpref,dpref,gpref,rc,rho0):
    return np.sqrt((disk(r,dpref))**2
               +(bulge(r,bpref))**2
               +(gas(r,gpref))**2
               +(halo(r,rc,rho0))**2)


def totalcurve(r,gpref,bpref,dpref,rc,rho0):
    return np.sqrt((gas(r,gpref)**2)
                   + (bulge(r,bpref)**2) 
                   + (disk(r,dpref)**2)
                   + (halo(r,rc,rho0)**2))


z=r_dat[0]
rr=[0,z]

# Define plotting function
def f(gpref,bpref,dpref,rc,rho0):
    
    # Define r
    xmax=20 #kpc
    rval = np.linspace(0,11.2,19)
    
    # Plot
    plt.figure(figsize=(11,6))
    plt.xlim(0,xmax)
    plt.ylim(0,360)
    
    plt.errorbar(r_dat,v_dat,yerr=v_err1,fmt='bo',label='Data')
    
    plt.plot(r_dat,bulge(r_dat,bpref),label=("Bulge"),color='orange')
    plt.plot(rr,[0,bulge(r_dat,bpref)[0]],color='orange') #straight lining connecting halo curve to origin, for visual aesthetic
    plt.plot(r_dat,disk(r_dat,dpref),label=("Disk"),color='purple')
    plt.plot(rr,[0,disk(r_dat,dpref)[0]],color='purple') #straight lining connecting halo curve to origin, for visual aesthetic
    plt.plot(r_dat,halo(r_dat,rc,rho0),label=("Halo"),color='green')
    plt.plot(rr,[0,halo(r_dat,rc,rho0)[0]],color='green') #straight lining connecting halo curve to origin, for visual aesthetic
    plt.plot(r_dat,gas(r_dat,gpref),label=("Gas"),color='blue')
    plt.plot(r_dat,totalcurve(r_dat,gpref,bpref,dpref,rc,rho0),label=("Total Curve"),color='red')
    plt.plot(rr,[0,totalcurve(r_dat,gpref,bpref,dpref,rc,rho0)[0]],color='red') #straight lining connecting halo curve to origin, for visual aesthetic

    plt.title(str("Interactive Rotation Curve - Galaxy: " + galaxy))
  
    plt.xlabel("Radius (kpc)")
    plt.ylabel("Velocity (km/s)")
    
    
    plt.legend(loc='lower right')
    
     
    #chi squared button commented out for now bc tiny graph
    residuals = v_dat - totalcurve(r_dat,gpref,bpref,dpref,rc,rho0)
    # Determining errors
    errors = v_err1**2 #second term is inclination uncertainty
    # Chi squared
    chisquared = np.sum(residuals**2/errors**2)
    #chisquared = stats.chisquare(v_dat,totalcurve(r,M,bpref,dpref,rc,rho0,gpref))
    reducedchisquared = chisquared * (1/(len(r_dat)-6))
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    plt.text(10,300,r"$\chi^2$: {:.5f}".format(chisquared)+'\n'+r"Reduced: {:.5f}".format(reducedchisquared),bbox=props)
    plt.annotate(str('Data source: '+ DataSource),xy=(0, 0), xytext=(0,5),xycoords=('axes fraction', 'figure fraction'),textcoords='offset points',size=10, ha='left', va='bottom')
    plt.show()

# Appearance
style = {'description_width': 'initial'}
layout = {'width':'600px'}

# Define slides
gpref = FloatSlider(min=0, max=5, step=0.1, 
                    value=gprefs,
                    description='Gas Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

bpref = FloatSlider(min=0, max=10, step=0.1, 
                    value=bprefs, 
                    description='Bulge Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

dpref = FloatSlider(min=0, max=5, step=0.1, 
                    value=dprefs, 
                    description='Disk Prefactor', 
                    readout_format='.2f', 
                    orientation='horizontal', 
                    style=style, layout=layout)

rc = FloatSlider(min=0, max=5, step=0.1, value=rcs, description='Halo Core Radius [kpc]', readout_format='.2f', orientation='horizontal', style=style, layout=layout)
#rc = fixed(best_rc)
#rho0 = fixed(best_rho0)

rho0 = FloatSlider(min=0, max=1e9, step=1e7, 
                    value=rho0s, 
                    description=r'Halo Surface Density [$M_{\odot} / pc^3$]', 
                    readout_format='.2e', 
                    orientation='horizontal', 
                    style=style, layout=layout)

# Interactive widget
def interactive_plot(f):
    interact = interactive(f,
                               gpref = gpref,
                               bpref = bpref, 
                               dpref = dpref, 
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
    bpref.value = bprefs
    dpref.value = dprefs
    rho0.value = rho0s
    rc.value = rcs
    gpref.value = gprefs

button.on_click(on_button_clicked)

# displaying button and its output together
VBox([button,out,interactive_plot(f)])

