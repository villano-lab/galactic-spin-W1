f=open('galaxy.txt', 'r')
f=f.read()

#for NGC5533:
Msun_B = 5.44         #constant, absolute magnitude of the sun in the B-band [source 4, Abstract]
Mabs_B = -21.22       #B-band magnitude of 5533 
Msun_Bsource= ' [unitless] (source 4, Abstract)'
Mabs_Bsource=' [unitless] (Source 2, Table 1)'
if f=='NGC5533':
    
    n = 2.7                #concentration parameter [unitless] (Source 2, Table A4)
    re=9.9                 #effective radius, [arcsec] (Source 2, table A4)
    D=54.3                 #distance, [Mpc] (Source 2, Table 1)
    L=5e9 
    q = 0.33               #intrinsic axis ratio [unitless] (Source 1, Table 1)
    #ML = 2.8              #mass-to-light ratio of bulge [unitless] (Source: 1, Table 1 for q = 0.33)
    i = 52 

    nsource=' [unitless] (Source 2, Table A4)'
    resource=' [arcsec] (Source 2, Table A4)'
    Dsource=' [Mpc] (Source 2, Table 1)'
    Lsource=' [Lsun] (Source 1, sec. 3, par. 1)'
    Msource=' [unitless] (Source 2, Table 1)'
    qsource=' [unitless] (Source 1, Table 1)'
    isource=' [degrees] (Source 2, Table A2)'
    
elif f=='NGC7814':
    
    #**************************
    #******NGC 7814 Bulge****
    #*************************

    # parameters for NGC 7814 provided in this paper: https://www.aanda.org/articles/aa/pdf/2011/07/aa16634-11.pdf 
    #parameters 
    n = 10       #concentration parameter [unitless] (Source 3, section 3.3 par.1)
    re = 2.16    #effective radius, [kpc] (Source 3, Table 4)
    L = 7e10     #bulge luminosity, [Lsun] (Source 3, Table 4)
    #ML = .71    #mass-to-light ratio of bulge given DM halo (Table 5)
    q = 0.61     #intrinsic axis ratio [unitless] (Source 3, Table 4)
    i=90         #inclination angle [degrees] (Source 3, Table 1)
    
    nsource=' [unitless] (Source 3, section 3.3 par.1)'
    resource=' [arcsec] (Source 3, Table 4)'
    Lsource=' [Lsun] (Source 3, Table 4)'
    Msource=' [unitless] (Source 2, Table 1)'
    qsource=' [unitless] (Source 3, Table 4)'
    isource=' [degrees] (Source 3, Table 1)'

elif f=='NGC891':
    
    #**************************
    #******NGC 891 Bulge****
    #*************************

    # parameters for NGC 891 provided in this paper: https://www.aanda.org/articles/aa/pdf/2011/07/aa16634-11.pdf 
    #parameters 
    n = 10        #concentration parameter (section 3.3 par.1)
    re = 1.8      #effective radius, [kpc] (table 4)
    L = 2.2e10    #bulge luminosity, [Lsun] (table 4)
    #ML = 1.63    #mass-to-light ratio of bulge given DM halo (table 5)
    D=9.5
    q = 0.68      #intrinsic axis ratio (table 4)
    i=89          #inclination angle [degrees] (table 1)

    nsource=' [unitless] (Source 3, section 3.3 par.1)'
    resource=' [arcsec] (Source 3, Table 4)'
    Lsource=' [Lsun] (Source 3, Table 4)'
    Msource=' [unitless] (Source 2, Table 1)'
    qsource=' [unitless] (Source 3, Table 4)'
    isource=' [degrees] (Source 3, Table 1)'
else:
    print("Oops! Make sure you typed your galaxy in correctly in cell 2!")