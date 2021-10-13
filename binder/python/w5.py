#!/usr/bin/env python
# coding: utf-8

# In[1]:


################################
########### Imports ############
################################
from sympy.utilities.lambdify import lambdify
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interactive, fixed, FloatSlider, HBox, Layout, Button, Label, Output, VBox
from IPython.display import display
import code7814 as wi
import matplotlib.image as mpimg
import sys
sys.path.append('../../python/')          # Defining path
import dataPython as dp
#getting so that the halo component/total RC curve updates as bh widgets update


def h(r,Mbh,arraysize):
    Me=wi.scale*Mbh*arraysize #units: [bh/dot]*[Msun/bh] * [dot] * [Msun] = Msun --> changes as sliders are moved
    MMM=Me*wi.viso #units: unitless * [Msun] = Msun
    return np.sqrt(MMM*wi.GG/wi.rcut) #is it rc or just r? the model looks better with rcut 
    #units: ( [Msun]* [(km/s)^2*(kpc/Msun)] / kpc )^.5 = km/s                                                                

def t(r,Mbh,arraysize):
    return np.sqrt((wi.d)**2
               +(wi.b)**2
               +(h(r,Mbh,arraysize))**2
               +(wi.g)**2)
img = mpimg.imread('../python/thing.jpg') #import special snowflake ngc 6814, which has visual diameter about 27.6kpc\ 

# Define plotting function
def f(arraysize,Mbh):
    
    #fractions: now changing the number of dots plottd within each bracket, which changes as the slider 
    #variable ("arraysize") changes
    
    #and trim to the first x elements of the pre-caclated radius arrays for each bracket
    arraysize=int(arraysize) #units: dot
    re = wi.rrr[:arraysize]
    anglee=wi.angle[:arraysize]
    
    #plot for each bracket, two equations for each bracket
    y=wi.c_y+re*np.sin(anglee) #y coordinates
    x=wi.c_x+re*np.cos(anglee) #x coordinates
    
    
    
    #below: changing the physical size of the plotted dots as the bh mass slider changes. the if else statement is fancy
    #but not necessary
    
    Mbh=1.5*Mbh #judging by eye
    if .1<=Mbh<1:
        Mbh=Mbh**1.5
    elif 1<=Mbh<2:
        Mbh=Mbh**2
    else:
        Mbh=Mbh**3 #display size of each marker
        
    
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
    f.set_figheight(15)
    f.set_figwidth(35)
    
     # Define r, for some reason this isn't working in code7814.py, so put directly here:
    data = dp.getXYdata_wXYerr('../testing/7814/ngc7814data')
    r = np.asarray(data['xx'])
    z=r[0]
    rr=[0,z]
    ax1.scatter(x,y,color="r",marker='o',s=Mbh)
    ax1.imshow(img)
    ax1.axis('off')
    ax2.tick_params(axis='both', which='major', labelsize=30)
    ax2.plot(r,h(r,Mbh,arraysize),label=("Dark Matter Halo"),color='green')
    ax2.plot(rr,[0,h(z,Mbh,arraysize)[0]],color='green') #straight lining connecting halo curve to origin, for visual aesthetic
    ax2.errorbar(r,wi.v_dat,yerr=wi.v_err1,fmt='bo',label='Data')
    ax2.plot(r,wi.b,label=("Bulge"),color='orange')
    ax2.plot(rr,[0,wi.b[0]],color='orange') #straight line for bulge
    ax2.plot(r,wi.d,label=("Disk"),color='purple')
    ax2.plot(rr,[0,wi.d[0]],color='purple') #straight line for disk
    ax2.plot(r,wi.g,label=("Gas"),color='blue')
    ax2.plot(r,t(r,Mbh,arraysize),label=("Total Curve"),color='red')
    ax2.plot(rr,[0,t(z,Mbh,arraysize)[0]],color='red') #straight line for total curve
    ax2.set_xlim([wi.minkpc, 19])
    ax2.set_ylim([0, 400])
    ax2.set_title('NGC 7814',fontsize = 40)
    ax2.set_ylabel('Velocity [km/s]',fontsize = 30)
    ax2.set_xlabel('radius [kpc]',fontsize = 30)
    ax2.legend(loc='upper right',prop={'size': 20},ncol=3)
    ax2.text(.1, 310, 'How close can you get? \nA reduced $\chi^2$ less than 1 is considered a good fit',fontsize=30)
    ax = plt.gca()
    xticks = ax.xaxis.get_major_ticks() 
    xticks[0].label1.set_visible(False)
    residuals = wi.v_dat - t(wi.r,Mbh,arraysize)
    # Determining errors
    errors = wi.v_err1**2 #second term is inclination uncertainty
    # Chi squared
    chisquared = np.sum(residuals**2/errors**2)
    #chisquared = stats.chisquare(v_dat,totalcurve(r,M,bpref,dpref,rc,rho00,gpref))
    reducedchisquared = chisquared * (1/(len(wi.r)-6))
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax2.text(0.4,360,r"Reduced $\chi^2$:"+'\n'+"{:.2f}".format(reducedchisquared),bbox=props,size=30)
    
   
    
style = {'description_width': 'initial'}
layout = {'width':'600px'}

################################
######## Define Sliders ########
################################

#mass of each black hole
Mbh = FloatSlider(min=wi.minmass, max=wi.maxmass, step=wi.minmass, 
                value=wi.start,
                description='Mass of each lil black holes [Msun]', 
                readout= True,
                readout_format='.1f', 
                orientation='horizontal', 
                style=style, layout=layout)

#number of projected black dots slider
arraysize = FloatSlider(min=wi.Min, max=wi.Max, step=wi.stepN, 
                value=wi.best_A, 
                description='Number of lil black holes multiplied by scale %.1e'%wi.scale, 
                readout= True,
                readout_format='.2d', 
                orientation='horizontal', 
                style=style, layout=layout)


def interactive_plot(f):
    interact = interactive(f, Mbh = Mbh, arraysize=arraysize,continuous_update=False)
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
    #display(Javascript('IPython.notebook.execute_cells_below()'))
    Mbh.value = wi.start
    arraysize.value=wi.best_A
button.on_click(on_button_clicked)

# displaying button and its output together
#VBox([button,out,interactive_plot(f)])


# In[ ]:




