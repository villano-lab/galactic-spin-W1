#!/usr/bin/env python
# coding: utf-8


Msun_B = 5.44          #constant, absolute magnitude of the sun in the B-band [source 4, Abstract]


n5533 = 2.7                #concentration parameter [unitless] (Source 2, Table A4)
re5533=9.9                 #effective radius, [arcsec] (Source 2, table A4)
D5533=54.3                 #distance, [Mpc] (Source 2, Table 1)
M_abs_B5533 = -21.22       #B-band magnitude of 5533 [unitless] (Source 2, Table 1)
q5533 = 0.33               #intrinsic axis ratio [unitless] (Source 1, Table 1)
#ML = 2.8              #mass-to-light ratio of bulge [unitless] (Source: 1, Table 1 for q = 0.33)
i5533 = 52 

nsource5533=' [unitless] (Source 2, Table A4)'
resource5533=' [arcsec] (Source 2, Table A4)'
Dsource5533=' [Mpc] (Source 2, Table 1)'
Msource5533=' [unitless] (Source 2, Table 1)'
qsource5533=' [unitless] (Source 1, Table 1)'
isource5533=' [degrees] (Source 2, Table A2)'




#**************************
#******NGC 7814 Bulge****
#*************************

# parameters for NGC 7814 provided in this paper: https://www.aanda.org/articles/aa/pdf/2011/07/aa16634-11.pdf 
#parameters 
n = 10       #concentration parameter (section 3.3 par.1)
re = 2.16    #effective radius, [kpc] (table 4)
L = 7e10     #bulge luminosity, [Lsun] (table 4)
ML = .71    #mass-to-light ratio of bulge given DM halo (table 5)
q = 0.61     #intrinsic axis ratio (table 4)
i=90         #inclination angle [degrees] (table 1)

#import GIPSY bulge for comparison bulge:
#bulge = np.loadtxt('../testing/7814reallybulge.dat')
prefactor=4.98 #from our fitting with NGC 7814
#radii,k,b = bulge.T #ignore k