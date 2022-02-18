### Imports ###
import numpy as np

# Read the file where the chosen galaxy name is located
f = open('python/galaxy_bulgeactivity.txt', 'r')
galaxy = f.read()

################
### NGC 5533 ###
################

if galaxy == 'NGC5533':
    check = True
    
    Mabs = -21.66     # B-band absolute magnitude of NGC 5533 (Norrdermeer, Van Der Hulst, 2007)
    n = 2.7           # Concentration parameter [unitless] (Noordermeer, Van Der Hulst, 2007)
    q = 0.33          # Intrinsic axis ratio [unitless] (Noordermeer, 2008)
    i = 52            # Inclination angle [unitless] (Noordermeer, Van Der Hulst, 2007)
    re_arcsec = 9.9   # Effective radius [arcsec] (Noordermeer, Van Der Hulst, 2007)
    Lb = 3.6e10       # Bulge luminosity [solar luminosity] (Calculated from absolute magnitude of the bulge)
    ML = 2.8          # Mass-to-light ratio of bulge [unitless] (Noordermeer, 2008)
    D_Mpc = 54.3      # Distance [Mpc] (Noordermeer, Van Der Hulst, 2007)

    Mabs_source =  ' [unitless] (Source #2, Table A4)'
    n_source =  ' [unitless] (Source #2, Table A4)'
    q_source =  ' [unitless] (Source #1, Table 1)'
    i_source =  ' [degrees] (Source #2, Table A2)'
    re_source = ' [kpc] (Source #2, Table A4)'
    Lb_source = ' [Lsun] (Source #2, calculated from absolute magnitude)'
    ML_source = ' [unitless] (Source #1, Table 1)'
    D_source =  ' [Mpc] (Source #2, Table 1)'
    
    # Unit conversion
    re_rad = re_arcsec * (np.pi / (3600 * 180))  # arcsec to radians
    D_kpc = D_Mpc * 1000                         # Mpc to kpc
    re_kpc = re_rad * D_kpc                      # Effective radius in kpc
    
###############
### NGC 891 ###
###############   

elif galaxy == 'NGC891':
    check = True
    
    # Source: https://www.aanda.org/articles/aa/pdf/2011/07/aa16634-11.pdf 
    Mabs = -21.12     # Absolute magnitude [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    n = 2.99          # Concentration parameter [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    q = 0.68          # Intrinsic axis ratio [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    i = 89            # Inclination angle [degrees] (Fraternali, Sancisi, and Kamphuis, 2011)
    re_kpc = 1.8      # Effective radius [kpc] (Fraternali, Sancisi, and Kamphuis, 2011)
    Lb = 2.2e10       # Bulge luminosity [solar luminosity] (Fraternali, Sancisi, and Kamphuis, 2011)
    ML = 1.63         # Mass-to-light ratio of bulge (Fraternali, Sancisi, and Kamphuis, 2011)
    D_Mpc = 9.5       # Distance [Mpc] (need this)
    
    Mabs_source =  ' [unitless] (Source #3, calculated from Luminosity)'
    n_source =  ' [unitless] (Source #3, Table 4)'
    q_source =  ' [unitless] (Source #3, Table 4)'
    i_source =  ' [degrees] (Source #3, Table 1)'
    re_source = ' [kpc] (Source #3, Table 4)'
    Lb_source = ' [Lsun] (Source #3, Table 4)'
    ML_source = ' [unitless] (Source #3, Table 5, with dark matter halo)'
    D_source =  ' [Mpc] (Source #3, Table 1)'
    
################
### NGC 7814 ###
################    

elif galaxy == 'NGC7814':
    check = True
    
    # Source: https://www.aanda.org/articles/aa/pdf/2011/07/aa16634-11.pdf 
    Mabs = -22.38     # Absolute magnitude [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    n = 4.0           # Concentration parameter [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    q = 0.61          # Intrinsic axis ratio [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    i = 90            # Inclination angle [degrees] (Fraternali, Sancisi, and Kamphuis, 2011)
    re_kpc = 2.16     # Effective radius [kpc] (Fraternali, Sancisi, and Kamphuis, 2011)
    Lb = 7e10         # Bulge luminosity [solar luminosity] (Fraternali, Sancisi, and Kamphuis, 2011)
    ML = 0.71         # Mass-to-light ratio of bulge (Fraternali, Sancisi, and Kamphuis, 2011)
    D_Mpc = 14.6      # Distance [Mpc] (need this)
    
    Mabs_source =  ' [unitless] (Source #3, calculated from Luminosity)'
    n_source =  ' [unitless] (Source #3, Table 4)'
    q_source =  ' [unitless] (Source #3, Table 4)'
    i_source =  ' [degrees] (Source #3, Table 1)'
    re_source = ' [kpc] (Source #3, Table 4)'
    Lb_source = ' [Lsun] (Source #3, Table 4)'
    ML_source = ' [unitless] (Source #3, Table 5, with dark matter halo)'
    D_source =  ' [Mpc] (Source #3, Table 1)'
    
else:
    check = False
    #print("Oops! Make sure you typed your galaxy in correctly in cell 2!")
    print("You have chosen the galaxy {}.".format(galaxy))
    
    
######################
### Check research ###
######################

from IPython.display import display, clear_output
import ipywidgets as widgets

button = widgets.Button(description="Reveal!")
def on_button_clicked(x):
    print(str('Galaxy '+ galaxy +':'))
    if check == True:
        print(str('Concentration parameter, n = ' + "{:.1f}".format(n) + n_source))
        print(str('Intrinsic axis ratio, q = ' + "{:.2f}".format(q) + q_source))
        print(str('Inclination angle, i = ' + "{:.1f}".format(i) + i_source))
        print(str('Effective radius, re = ' + "{:.1f} ".format(re_kpc) + re_source))
        print(str('Luminosity of the bulge, L = ' + "{:.1e}".format(Lb) + Lb_source))
        print(str('Distance, D = ' + "{:.1f}".format(D_Mpc) + D_source))
    if check == False:
        print("Sorry, the parameters for your chosen galaxy is not in our library.")
button.on_click(on_button_clicked)