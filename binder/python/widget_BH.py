#########################
### Black Hole Widget ###
#########################

###############
### Imports ###
###############

import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interactive, fixed, FloatSlider, HBox, Layout, Button, Label, Output, VBox
from IPython.display import display
import components as comp
import dataPython as dp

####################
### Galaxy image ###
####################

# Import galaxy image
img = plt.imread("images/A_spiral_snowflake.jpg")            # Import image of NGC 6814

# Find the center by eye from the image
c_x = 1960
c_y = 1800

# Kpc limits, visual guess based on the galaxy in the image chosen:
minkpc = 0
maxkpc = 100

##################
### Parameters ###
##################

# units: scale = [#number of actual black holes / plotted dot]
kpctopixels = 20                # visual scaling, varies depending on size of galaxy image (and actual size of galaxy)
r1 = minkpc*kpctopixels
r2 = maxkpc*kpctopixels

# For mass of black hole slider:
minmassBH = 0.1                 # solar masses, arbitrary
maxmassBH = 3.8                 # solar masses, size of the smallest black hole ever discovered according to
                                # https://www.scientificamerican.com/gallery/the-smallest-known-black-hole/
defaultmass = 1.5               # default mass value for slider
stepM = minmassBH               # step of mass of bh's

### NGC 5533 ### 
# For number of black holes slider
scale_5533 = 1e6                      # scale is neccessary to be a constant, otherwise the widget will freeze up the computer! 
                                      # scale is how many black holes represent each dot on image
minnumberBH_5533 = 50                 # min number of black holes (this is multiplied by the scale)
maxnumberBH_5533 = 5e2                # max number of black holes (this is multiplied by the scale)
defaultnumber_5533 = 245              # default number of bh's for slider (this is multiplied by the scale)
stepN_5533 = 5                        # step of # of bh's for slider (this is multiplied by the scale)

# For cutoff radius
minrcutBH_5533 = 0.1                  # min cutoff radius for black holes
maxrcutBH_5533 = 3.0                  # max cutoff radius for black holes
defaultrcut_5533 = 1.4                # default cutoff radius for black holes
stepRCUT_5533 = minrcutBH_5533        # step of cutoff radius for black holes

# Generate random positions for black holes inside the donut
rand_radius_5533 = lambda rad1,rad2: np.random.uniform(rad1,rad2,int(maxnumberBH_5533))
rand_angle_5533 = np.random.uniform(0,2*np.pi,int(maxnumberBH_5533))  # angle 0 to 360 degrees for full circle (donut) for each bracket


### NGC 7814 ###
# For number of black holes slider
scale_7814 = 1e6                      # scale is neccessary to be a constant, otherwise the widget will freeze up the computer! 
                                      # scale is how many black holes represent each dot on image
minnumberBH_7814 = 1                  # min number of black holes (this is multiplied by the scale)
maxnumberBH_7814 = 1000               # max number of black holes (this is multiplied by the scale)
defaultnumber_7814 = 600              # default number of bh's for slider (this is multiplied by the scale)
stepN_7814 = 5                        # step of # of bh's for slider (this is multiplied by the scale)

# For cutoff radius
minrcutBH_7814 = 0.1                  # min cutoff radius for black holes
maxrcutBH_7814 = 4.0                  # max cutoff radius for black holes
defaultrcut_7814 = 1.6                # default cutoff radius for black holes
stepRCUT_7814 = minrcutBH_7814        # step of cutoff radius for black holes

# Generate random positions for black holes inside the donut
rand_radius_7814 = lambda rad1,rad2: np.random.uniform(rad1,rad2,int(maxnumberBH_7814))
rand_angle_7814 = np.random.uniform(0,2*np.pi,int(maxnumberBH_7814))  # angle 0 to 360 degrees for full circle (donut) for each bracket

style = {'description_width': 'initial'}
layout = {'width':'800px'}

####################################
### Plotting function for widget ###
####################################

### NGC 5533 ###
def f5533(arraysize,massMiniBH,rcut):
    
    _,fit_dict = comp.bestfit(comp.totalvelocity_miniBH,'NGC5533')
    
    bpref = fit_dict['bpref']
    dpref = fit_dict['dpref']
    gpref = fit_dict['gpref']
    Mbh = fit_dict['Mbh']
    
    # Define radius
    r = np.linspace(
        np.min(comp.galdict('NGC5533')['m_radii']), 
        np.max(comp.galdict('NGC5533')['m_radii']),500)
    
    # Change input to an integer
    arraysize = int(arraysize)     # units: dot
    
    # Random radii for the galaxy image
    # Trim to the first x elements (arraysize) of the pre-calculated radius arrays for each bracket
    radius_trim = rand_radius_5533(r1,r2)[:arraysize]
    angle_trim = rand_angle_5533[:arraysize]
    
    # x and y coordinates for plotting
    x = c_x + radius_trim*np.cos(angle_trim)     # x coordinates
    y = c_y + radius_trim*np.sin(angle_trim)     # y coordinates
    
    # Random radii for the plot
    radius_plot = rand_radius_5533(minkpc,maxkpc)[:arraysize]
    radius_plot = np.sort(radius_plot)           # sort array
    radius_plot[0] = 0.2    # first element of the radius should be close to zero because spline function won't plot up to zero
      
    # Set up two plots next to each other
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
    plt.subplots_adjust(wspace=0, hspace=0)
    f.set_figheight(12)
    f.set_figwidth(32)
    
    # Changing the size of each dot as the mass of the black hole changes
    # Without this the dots would be either too small or too big  
    dotsize = np.linspace(5,12,int(maxmassBH/stepM + 1))            # array of sizes from 5 to 12 for dots in scatterplot
    mass = np.round(np.arange(minmassBH,maxmassBH+stepM,stepM),1)   # array of masses, rounded to 1 decimal
    loc = np.where(mass == massMiniBH)[0]                           # which element equals to the mBH value, returns position
    BHsize = dotsize[loc[0]]                                        # picks out a dotsize for that mass
    
    # First plot - image with black holes
    ax1.plot(x,y, linestyle='none', markerfacecolor='orangered', marker="o", markeredgecolor="maroon", markersize=BHsize)
    ax1.imshow(img)
    ax1.set_title("Each dot representing 1 million tiny black holes.", fontsize=25)
    ax1.set_xlim(0,3970)
    ax1.set_ylim(0,3970)
    ax1.axis('off')
    
    # Second plot - rotation curve   
    ax2.plot(r,comp.halo_BH(r,scale_5533,arraysize,massMiniBH,rcut),
             label=("Dark Matter Halo - Tiny Black Holes"),color='green')
    ax2.errorbar(comp.galdict('NGC5533')['m_radii'],comp.galdict('NGC5533')['m_velocities'],
                 yerr=comp.galdict('NGC5533')['m_v_errors'],fmt='bo',label='Data')
    ax2.plot(r,comp.blackhole(r,Mbh,'NGC5533'),label=("Central Black Hole"),color='black')
    ax2.plot(r,comp.bulge(r,bpref,'NGC5533'),label=("Bulge"),color='orange')
    ax2.plot(r,comp.disk(r,dpref,'NGC5533'),label=("Disk"),color='purple')
    ax2.plot(r,comp.gas(r,gpref,'NGC5533'),label=("Gas"),color='blue')
    ax2.plot(r,comp.totalvelocity_miniBH(r,scale_5533,arraysize,massMiniBH,rcut,
                                  bpref,dpref,gpref,
                                  Mbh,'NGC5533'),label=("Total Curve"),color='red')
    ax2.set_title('NGC 5533',fontsize=40)
    ax2.set_ylabel('Velocity [km/s]',fontsize=25)
    ax2.set_xlabel('Radius [kpc]',fontsize=25)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.set_xlim(0,100)
    ax2.set_ylim(0,400)
    ax2.legend(bbox_to_anchor=(1,1), loc="upper left", fontsize=20) 

    # Residuals
    residuals = comp.galdict('NGC5533')['m_velocities'] 
    - comp.totalvelocity_miniBH(comp.galdict('NGC5533')['m_radii'],
                                scale_5533,arraysize,massMiniBH,rcut,
                                                        bpref, dpref, gpref, Mbh,
                                                        'NGC5533')
    # Chi squared
    chisquared = np.sum(residuals**2/comp.galdict('NGC5533')['m_v_errors']**2)
    dof = len(comp.galdict('NGC5533')['m_radii']) - 6       # number of degrees of freedom = number of observed data - number of fitting parameters
    reducedchisquared = chisquared / dof
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax2.text(98,390,r"Reduced $\chi^2$: {:.2f}".format(reducedchisquared),ha='right',va='top',bbox=props,size=22)
        
### NGC 7814 ###
def f7814(arraysize,massMiniBH,rcut):
    
    _,fit_dict = comp.bestfit(comp.totalvelocity_miniBH,'NGC7814')
    
    bpref = fit_dict['bpref']
    dpref = fit_dict['dpref']
    gpref = fit_dict['gpref']
    
    # Define radius
    r = np.linspace(
        np.min(comp.galdict('NGC7814')['m_radii']),
        np.max(comp.galdict('NGC7814')['m_radii']),500)
    
    # Change input to an integer
    arraysize = int(arraysize)     # units: dot
    
    # Random radii for the galaxy image
    # Trim to the first x elements (arraysize) of the pre-calculated radius arrays for each bracket
    radius_trim = rand_radius_7814(r1,r2)[:arraysize]
    angle_trim = rand_angle_7814[:arraysize]
    
    # x and y coordinates for plotting
    x = c_x + radius_trim*np.cos(angle_trim)     # x coordinates
    y = c_y + radius_trim*np.sin(angle_trim)     # y coordinates
    
    # Random radii for the plot
    radius_plot = rand_radius_7814(minkpc,maxkpc)[:arraysize]
    radius_plot = np.sort(radius_plot)           # sort array
    radius_plot[0] = 0.2    # first element of the radius should be close to zero because spline function won't plot up to zero
      
    # Set up two plots next to each other
    f, (ax3, ax4) = plt.subplots(1, 2, sharey=False)
    plt.subplots_adjust(wspace=0, hspace=0)
    f.set_figheight(12)
    f.set_figwidth(32)
    
    # Changing the size of each dot as the mass of the black hole changes
    # Without this the dots would be either too small or too big  
    dotsize = np.linspace(5,12,int(maxmassBH/stepM + 1))            # array of sizes from 5 to 12 for dots in scatterplot
    mass = np.round(np.arange(minmassBH,maxmassBH+stepM,stepM),1)   # array of masses, rounded to 1 decimal
    loc = np.where(mass == massMiniBH)[0]                           # which element equals to the mBH value, returns position
    BHsize = dotsize[loc[0]]                                        # picks out a dotsize for that mass
    
    # First plot - image with black holes
    ax3.plot(x,y, linestyle='none', markerfacecolor='orangered', marker="o", 
             markeredgecolor="maroon", markersize=BHsize)
    ax3.imshow(img)
    ax3.set_title("Each dot representing 1 million tiny black holes.", fontsize=25)
    ax3.set_xlim(0,3970)
    ax3.set_ylim(0,3970)
    ax3.axis('off')
    
    # Second plot - rotation curve   
    ax4.plot(r,comp.halo_BH(r,scale_7814,arraysize,massMiniBH,rcut),
             label=("Dark Matter Halo - Tiny Black Holes"),color='green')
    ax4.errorbar(comp.galdict('NGC7814')['m_radii'],
                 comp.galdict('NGC7814')['m_velocities'],yerr=comp.galdict('NGC7814')['m_v_errors'],fmt='bo',label='Data')
    ax4.plot(r,comp.bulge(r,bpref,'NGC7814'),label=("Bulge"),color='orange')
    ax4.plot(r,comp.disk(r,dpref,'NGC7814'),label=("Disk"),color='purple')
    ax4.plot(r,comp.gas(r,gpref,'NGC7814'),label=("Gas"),color='blue')
    ax4.plot(r,comp.totalvelocity_miniBH(r,scale_7814,arraysize,massMiniBH,
                                                  rcut,bpref,dpref,
                                                  gpref,0,'NGC7814'),label=("Total Curve"),color='red')
    ax4.set_title('NGC 7814',fontsize=40)
    ax4.set_ylabel('Velocity [km/s]',fontsize=25)
    ax4.set_xlabel('Radius [kpc]',fontsize=25)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)
    ax4.set_xlim(0,np.max(comp.galdict('NGC7814')['m_radii']))
    ax4.set_ylim(0,400)
    ax4.legend(bbox_to_anchor=(1,1), loc="upper left", fontsize=20) 

    # Residuals
    residuals = comp.galdict('NGC7814')['m_velocities'] - comp.totalvelocity_miniBH(
        comp.galdict('NGC7814')['m_radii'], 
        scale_7814,arraysize,massMiniBH,rcut, bpref, dpref, gpref, 0, 'NGC7814')
    # Chi squared
    chisquared = np.sum(residuals**2/comp.galdict('NGC7814')['m_v_errors']**2)
    dof = len(comp.galdict('NGC5533')['m_radii']) - 5       # number of degrees of freedom = number of observed data - number of fitting parameters
    reducedchisquared = chisquared / dof  
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax4.text(19,390,r"Reduced $\chi^2$: {:.2f}".format(reducedchisquared),ha='right',va='top',bbox=props,size=22)
    
################################
######## Define Sliders ########
################################

### NGC 5533 ###
# Number of projected black dots slider
arraysize_5533 = FloatSlider(min=minnumberBH_5533, max=maxnumberBH_5533, step=stepN_5533, 
                value=defaultnumber_5533, 
                description='Number of millions of tiny black holes (increasing by {:.0f} million)'.format(stepN_5533), 
                readout=True,
                readout_format='.2d', 
                orientation='horizontal', 
                style=style, layout=layout)

# Mass of each black hole
massMiniBH_5533 = FloatSlider(min=minmassBH, max=maxmassBH, step=stepM, 
                value=defaultmass,
                description='Mass of each tiny black hole (in solar masses, increasing by {:.1f})'.format(stepM), 
                readout=True,
                readout_format='.1f',
                orientation='horizontal', 
                style=style, layout=layout)

# Cutoff radius
rcut_5533 = FloatSlider(min=minrcutBH_5533, max=maxrcutBH_5533, step=stepRCUT_5533, 
                value=defaultrcut_5533,
                description='Cutoff radius (in kpc, increasing by {:.1f})'.format(stepRCUT_5533), 
                readout=True,
                readout_format='.1f',
                orientation='horizontal', 
                style=style, layout=layout)

### NGC 7814 ###
# Number of projected black dots slider
arraysize_7814 = FloatSlider(min=minnumberBH_7814, max=maxnumberBH_7814, step=stepN_7814, 
                value=defaultnumber_7814, 
                description='Number of millions of tiny black holes (increasing by {:.0f} million)'.format(stepN_7814), 
                readout=True,
                readout_format='.2d', 
                orientation='horizontal', 
                style=style, layout=layout)

# Mass of each black hole
massMiniBH_7814 = FloatSlider(min=minmassBH, max=maxmassBH, step=stepM, 
                value=0.5,
                description='Mass of each tiny black hole (in solar masses, increasing by {:.1f})'.format(stepM), 
                readout=True,
                readout_format='.1f',
                orientation='horizontal', 
                style=style, layout=layout)

# Cutoff radius
rcut_7814 = FloatSlider(min=minrcutBH_7814, max=maxrcutBH_7814, step=stepRCUT_7814, 
                value=defaultrcut_7814,
                description='Cutoff radius (in kpc, increasing by {:.1f})'.format(stepRCUT_7814), 
                readout=True,
                readout_format='.1f',
                orientation='horizontal', 
                style=style, layout=layout)


def interactive_plot_5533(f5533):
    interact = interactive(f5533, arraysize=arraysize_5533, 
                           scale=scale_5533, massMiniBH=massMiniBH_5533, 
                           rcut=rcut_5533, continuous_update=False)
    return interact

def interactive_plot_7814(f7814):
    interact = interactive(f7814, arraysize=arraysize_7814, 
                           scale=scale_7814, massMiniBH=massMiniBH_7814, 
                           rcut=rcut_7814, continuous_update=False)
    return interact

################################
########### Button #############
################################

### NGC 5533 ###
# Button to revert back to Best Fit
button_5533 = Button(
    description="Best Fit",
    button_style='warning', # 'success', 'info', 'warning', 'danger' or ''
    icon='check')
out_5533 = Output()

def on_button_clicked_5533(_):
    arraysize_5533.value = defaultnumber_5533
    massMiniBH_5533.value = defaultmass
    rcut_5533.value = defaultrcut_5533
button_5533.on_click(on_button_clicked_5533)

### NGC 7814 ###
# Button to revert back to Best Fit
button_7814 = Button(
    description="Best Fit",
    button_style='warning', # 'success', 'info', 'warning', 'danger' or ''
    icon='check')
out_7814 = Output()

def on_button_clicked_7814(_):
    arraysize_7814.value = defaultnumber_7814
    massMiniBH_7814.value = 0.5
    rcut_7814.value = defaultrcut_7814
button_7814.on_click(on_button_clicked_7814)