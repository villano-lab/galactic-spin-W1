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
import NGC5533_components as comp5533
import ONE_widget_library as comp7814
import dataPython as dp

####################
### Galaxy image ###
####################

# Import galaxy image
img = plt.imread("images/A_spiral_snowflake.jpg")            # import special snowflake ngc 6814, 
                                                             # which has visual diameter about 27.6 kpc    
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
kpctopixels = 20       # visual scaling, varies depending on size of galaxy image (and actual size of galaxy)
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
defaultnumber_5533 = 2.2e2            # default number of bh's for slider (this is multiplied by the scale)
stepN_5533 = 10                       # step of # of bh's for slider (this is multiplied by the scale)

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
scale_7814 = 2.5e7                    # scale is neccessary to be a constant, otherwise the widget will freeze up the computer! 
                                      # scale is how many black holes represent each dot on image
minnumberBH_7814 = 1                  # min number of black holes (this is multiplied by the scale)
maxnumberBH_7814 = 100                # max number of black holes (this is multiplied by the scale)
defaultnumber_7814 = 50               # default number of bh's for slider (this is multiplied by the scale)
stepN_7814 = 10                       # step of # of bh's for slider (this is multiplied by the scale)

# For cutoff radius
minrcutBH_7814 = 0.1                  # min cutoff radius for black holes
maxrcutBH_7814 = 4.0                  # max cutoff radius for black holes
defaultrcut_7814 = 2.1                # default cutoff radius for black holes
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
    f.set_figheight(10)
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
    ax1.set_title("Each dot representing 1 million tiny black holes.", fontsize=20)
    ax1.set_xlim(0,3970)
    ax1.set_ylim(0,3970)
    ax1.axis('off')
    
    # Second plot - rotation curve   
    ax2.plot(radius_plot,comp5533.halo_BH(radius_plot,scale_5533,arraysize,massMiniBH,rcut),
             label=("Dark Matter Halo - Tiny Black Holes"),color='green')
    ax2.errorbar(comp5533.r_dat,comp5533.v_dat,yerr=comp5533.v_err1,fmt='bo',label='Data')
    ax2.plot(radius_plot,comp5533.blackhole(radius_plot,comp5533.best_Mbh),label=("Central Black Hole"),color='black')
    ax2.plot(radius_plot,comp5533.bulge(radius_plot,comp5533.best_bpref),label=("Bulge"),color='orange')
    ax2.plot(radius_plot,comp5533.disk(radius_plot,comp5533.best_dpref),label=("Disk"),color='purple')
    ax2.plot(radius_plot,comp5533.gas(radius_plot,comp5533.best_gpref),label=("Gas"),color='blue')
    ax2.plot(radius_plot,comp5533.totalvelocity(radius_plot,scale_5533,arraysize,massMiniBH,rcut,
                                                comp5533.best_Mbh,comp5533.best_bpref,comp5533.best_dpref,comp5533.best_gpref),
                                                 label=("Total Curve"),color='red')
    ax2.set_title('NGC 5533',fontsize = 40)
    ax2.set_ylabel('Velocity [km/s]',fontsize = 25)
    ax2.set_xlabel('radius [kpc]',fontsize = 25)
    ax2.tick_params(axis='x', labelsize=16)
    ax2.tick_params(axis='y', labelsize=16)
    ax2.set_xlim(0,100)
    ax2.set_ylim(0,400)
    
    ax2.legend(bbox_to_anchor=(1,1), loc="upper left") 
    #ax2.legend(loc='upper center',prop={'size': 16},ncol=3)

    # Residuals
    residuals = comp5533.v_dat - comp5533.totalvelocity(comp5533.r_dat,scale_5533,arraysize,massMiniBH,rcut,
                                                        comp5533.best_Mbh,comp5533.best_bpref,
                                                        comp5533.best_dpref,comp5533.best_gpref)
    # Chi squared
    chisquared = np.sum(residuals**2/comp5533.v_err1**2)
    dof = len(comp5533.r_dat) - 6       # number of degrees of freedom = number of observed data - number of fitting parameters
    reducedchisquared = chisquared / dof
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax2.text(73,330,r"Reduced $\chi^2$:"+'\n'+"{:.2f}".format(reducedchisquared),ha='left',va='top',bbox=props,size=20)
       
        
### NGC 7814 ###
def f7814(arraysize,massMiniBH,rcut):
    
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
    f.set_figheight(10)
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
    ax3.set_title("Each dot representing 1 million tiny black holes.", fontsize=20)
    ax3.set_xlim(0,3970)
    ax3.set_ylim(0,3970)
    ax3.axis('off')
    
    # Second plot - rotation curve   
    ax4.plot(radius_plot,comp7814.halo_BH(radius_plot,scale_7814,arraysize,massMiniBH,rcut),
             label=("Dark Matter Halo - Tiny Black Holes"),color='green')
    ax4.errorbar(comp7814.r_dat,comp7814.v_dat,yerr=comp7814.v_err1,fmt='bo',label='Data')
    ax4.plot(radius_plot,comp7814.bulge(radius_plot,comp7814.bprefs),label=("Bulge"),color='orange')
    ax4.plot(radius_plot,comp7814.disk(radius_plot,comp7814.dprefs),label=("Disk"),color='purple')
    ax4.plot(radius_plot,comp7814.gas(radius_plot,comp7814.gprefs),label=("Gas"),color='blue')
    ax4.plot(radius_plot,comp7814.totalvelocityBH(radius_plot,scale_7814,arraysize,massMiniBH,
                                                  rcut,comp7814.bprefs,comp7814.dprefs,
                                                  comp7814.gprefs),label=("Total Curve"),color='red')
    ax4.set_title('NGC 7814',fontsize = 40)
    ax4.set_ylabel('Velocity [km/s]',fontsize = 25)
    ax4.set_xlabel('radius [kpc]',fontsize = 25)
    ax4.tick_params(axis='x', labelsize=16)
    ax4.tick_params(axis='y', labelsize=16)
    ax4.set_xlim(0,np.max(comp7814.r_dat))
    ax4.set_ylim(0,400)
    ax4.legend(bbox_to_anchor=(1,1), loc="upper left") 

    # Residuals
    residuals = comp7814.v_dat - comp7814.totalvelocityBH(comp7814.r_dat,scale_7814,arraysize,massMiniBH,rcut,
                                                          comp7814.bprefs,comp7814.dprefs,comp7814.gprefs)
    # Chi squared
    chisquared = np.sum(residuals**2/comp7814.v_err1**2)
    dof = len(comp5533.r_dat) - 4       # number of degrees of freedom = number of observed data - number of fitting parameters
    reducedchisquared = chisquared / dof  
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax4.text(14.5,350,r"Reduced $\chi^2$:"+'\n'+"{:.2f}".format(reducedchisquared),ha='left',va='top',bbox=props,size=20)

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
                value=defaultmass,
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
    arraysize_5533.value = defaultnumber
    massMiniBH_5533.value = defaultmass
    rcut_5533.value = defaultrcut
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
    massMiniBH_7814.value = defaultmass_7814
    rcut_7814.value = defaultrcut_7814
button_7814.on_click(on_button_clicked_7814)