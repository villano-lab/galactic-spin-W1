"""A collection of dictionaries containing data for various galaxies.

This module exists to support other modules by organizing data in a quasi-standardized manner.
"""
#Data to fit to for each galaxy to be used in workshop

###############################
########## Imports ############
###############################

import dataPython         as dp
import numpy              as np
import scipy.interpolate  as inter

###############################
########### NGC5533 ###########
###############################

NGC5533 = {
    # Load data from files for Noordermeer 2008 band and fitted curves
    
    # 'raw' in the sense that I haven't split everything up. Idk there's probably a better way to name these
    'raw_total'        : dp.getXYdata('data/NGC5533/noord-120kpc-total.txt'     ),
    'raw_blackhole'    : dp.getXYdata('data/NGC5533/noord-120kpc-blackhole.txt' ),
    'raw_bulge'        : dp.getXYdata('data/NGC5533/noord-120kpc-bulge.txt'     ),
    'raw_disk'         : dp.getXYdata('data/NGC5533/noord-120kpc-disk.txt'      ),
    'raw_halo'         : dp.getXYdata('data/NGC5533/noord-120kpc-halo.txt'      ),
    'raw_gas'          : dp.getXYdata('data/NGC5533/noord-120kpc-gas.txt'       ),
    'raw_band_btm'     : dp.getXYdata('data/NGC5533/noord-120kpc-bottomband.txt'),
    'raw_band_top'     : dp.getXYdata('data/NGC5533/noord-120kpc-topband.txt'   ),
    
    # Get data from 100kpc file
    'measured_data'    : dp.getXYdata_wXYerr('data/NGC5533/100kpc_data.txt'),
    
    # Some constants
    'i'                : 52,    # Inclination angle [degrees] (Fraternali, Sancisi, and Kamphuis, 2011)  
    'D_Mpc'            : 54.3,  # Distance [Mpc]
    'Mabs'             : -21.66,# B-band absolute magnitude of NGC 5533 (Norrdermeer, Van Der Hulst, 2007)
    'check'            : True
}
"""Data for galaxy NGC5533. Parameters and data measurements in this dictionary are from [Noordermeer2007]_ unless noted otherwise.

:keys:
    blackhole: [dict] Further information pertaining to the black hole component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the black hole component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the black hole component.
        
        Mbh: [float] Mass of the black hole.
        
        r: [array] Radial (kpc) values of the black hole component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the black hole component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the black hole component.

        v: [array] Velocity (km/s) values of the black hole component.
    
    bulge: [dict] Further information pertaining to the bulge component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the bulge component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the bulge component.

        Lb: [float] Luminosity (solar luminosities) of the bulge.

        ML: [float] Mass-to-light ratio of the bulge.

        n: [float] Concentration parameter.

        q: [float] Intrinsic axis ratio.

        r: [array] Radial (kpc) values of the bulge component.

        re_arcsec: [float] Effective radius (arcsec).

        re_rad: [float] Effective radius (radians).

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the bulge component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the bulge component.

        v: [array] Velocity (km/s) values of the bulge component.

    check: [bool] A flag indicating that the necessary data for the "Literature Search Check" in `07_Bonus_Bulge_Rotation_Curve.ipynb <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/07_Bonus_Bulge_Rotation_Curve.ipynb>`_ is available.

    disk: [dict] Further information pertaining to the disk component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the disk component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the disk component.
        
        r: [array] Radial (kpc) values of the disk component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the disk component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the disk component.

        v: [array] Velocity (km/s) values of the disk component.

    D_Mpc: [float] The distance to the galaxy in Megaparsecs.

    D_kpc: [float] The distance to the galaxy in kiloparsecs.

    galaxyname: [string] Name of the galaxy.

    gas: [dict] Further information pertaining to the gas component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the gas component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the gas component.
        
        r: [array] Radial (kpc) values of the gas component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the gas component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the gas component.

        v: [array] Velocity (km/s) values of the gas component.

    halo: [dict] Further information pertaining to the halo component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the halo component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the halo component.
        
        r: [array] Radial (kpc) values of the halo component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the halo component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the halo component.

        v: [array] Velocity (km/s) values of the halo component.

    i: [int] The inclination angle, in degrees, of the galaxy.

    m_r_errors: [array] Errors of measured radii (kpc).

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).

    Mabs: [float] The absolute magnitude of the galaxy.

    massBH: [float] Mass of central black hole (solar masses).

    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ex: [list] Float values representing the error on the radial (kpc) data.

        ey: [list] Float values representing the error on the velocity (km/s) data.

    n_band_btm: [scipy.interpolate._bsplines.BSpline] Spline of the bottom side of the band.

    n_band_top: [scipy.interpolate._bsplines.BSpline] Spline of the top side of the band.

    n_cb: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the bottom side of the error band.

    n_ct: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the top side of the error band.

    n_kb: [int] Degree of spline returned by using scipy.interpolate.splrep on the bottom side of the error band.

    n_kt: [int] Degree of spline returned by using scipy.interpolate.splrep on the top side of the error band.

    n_r_btmband: [array] Radial values of the bottom side of the error band.

    n_r_topband: [array] Radial values of the top side of the error band.

    n_tb: [array] Vector of knots returned by using scipy.interpolate.splrep on the bottom side of the error band.

    n_tt: [array] Vector of knots returned by using scipy.interpolate.splrep on the top side of the error band.

    n_v_bandwidth: [array] Width of the error band, in km/s.

    n_v_btmband: [array] Velocity values of the bottom side of the error band.

    n_v_btmband: [array] Velocity values of the top side of the error band.

    raw_band_btm: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    raw_band_top: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.
    
    raw_blackhole: [dict] Data representing the black hole's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    raw_bulge: [dict] Data representing the bulge's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    raw_disk: [dict] Data representing the disk's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.
    
    raw_gas: [dict] Data representing the gas' contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    raw_halo: [dict] Data representing the halo's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.    
    
    raw_total: [dict] Data representing the total theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    rc: [float] Core radius (kpc).

    rho0: [float] Central mass density (solar mass/kpc^3).

    sources: [dict] Printing helper variables indicating units and sources for certain variables.
        D: [string] Units and source for distance to galaxy.
        
        i: [string] Units and source for inclination angle.

        Lb: [string] Units and source for bulge luminosity.
        
        Mabs: [string] Units and source for absolute magnitude.

        ML: [string] Units and source for bulge mass-light ratio.

        n: [string] Units and source for concentration parameter.

        re: [string] Units and source for effective radius.

        q: [string] Units and source for intrinsic axis ratio.

    total: [dict] Further information pertaining to the total theoretical rotation curve.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the total curve.
        
        k: [int] Degree of the spline returned by using scipy.interpolate.splrep on the total curve.

        r: [array] Radial (kpc) values of the total curve.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the total curve.

        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the total curve.

        v: [array] Velocity (km/s) values of the total curve.

"""

# Parameters ########################
NGC5533['galaxyname'] = 'NGC 5533'   # NGC catalog number of the galaxy
NGC5533['rho0'] = 0.31e9       # central mass density (in solar mass/kpc^3), Source: Noordermeer (2007)     
NGC5533['rc'] = 1.4            # core radius (in kpc), Source: Noordermeer (2007)
NGC5533['massbh'] = 2.7e9      # mass of central black hole (in solar masses), Source: Noordermeer (2007)

#Organize 100kpc data
NGC5533['m_radii']      = np.asarray(NGC5533['measured_data']['xx'])
NGC5533['m_velocities'] = np.asarray(NGC5533['measured_data']['yy'])
NGC5533['m_r_errors']   = np.asarray(NGC5533['measured_data']['ex'])
NGC5533['m_v_errors']   = np.asarray(NGC5533['measured_data']['ey'])

#Organize band data
NGC5533['n_r_btmband']   = np.asarray(NGC5533['raw_band_btm']['xx'])
NGC5533['n_v_btmband']   = np.asarray(NGC5533['raw_band_btm']['yy'])
NGC5533['n_r_topband']   = np.asarray(NGC5533['raw_band_top']['xx'])
NGC5533['n_v_topband']   = np.asarray(NGC5533['raw_band_top']['yy'])
NGC5533['n_v_bandwidth'] = (NGC5533['n_v_topband'] - NGC5533['n_v_btmband'])/2
NGC5533['n_v_bandwidth'] = NGC5533['n_v_bandwidth'][0::28] #For weights, v_errors and band must line up.
NGC5533['n_v_bandwidth'] = NGC5533['n_v_bandwidth'][1:]
    
# Smoothing
NGC5533['n_tb'], NGC5533['n_cb'], NGC5533['n_kb'] = inter.splrep(NGC5533['n_r_btmband'],NGC5533['n_v_btmband'])
NGC5533['n_tt'], NGC5533['n_ct'], NGC5533['n_kt'] = inter.splrep(NGC5533['n_r_topband'],NGC5533['n_v_topband'])
NGC5533['n_band_btm'] = inter.BSpline(NGC5533['n_tb'], NGC5533['n_cb'], NGC5533['n_kb'])
NGC5533['n_band_top'] = inter.BSpline(NGC5533['n_tt'], NGC5533['n_ct'], NGC5533['n_kt'])

# Total Curve #######################
NGC5533['total'] = {
    'r' : np.asarray(NGC5533['raw_total']['xx']),
    'v' : np.asarray(NGC5533['raw_total']['yy'])
}
NGC5533['total']['t'], NGC5533['total']['c'], NGC5533['total']['k'] = inter.splrep(NGC5533['total']['r'], NGC5533['total']['v'])
NGC5533['total']['spline'] = inter.BSpline(NGC5533['total']['t'], NGC5533['total']['c'], NGC5533['total']['k'])

# Black Hole ########################
NGC5533['blackhole'] = {
    'r'  : np.asarray(NGC5533['raw_blackhole']['xx']),
    'v'  : np.asarray(NGC5533['raw_blackhole']['yy']),
    'Mbh': 2.7e9
}
NGC5533['blackhole']['t'], NGC5533['blackhole']['c'], NGC5533['blackhole']['k'] = inter.splrep(NGC5533['blackhole']['r'], NGC5533['blackhole']['v'])
NGC5533['blackhole']['spline'] = inter.BSpline(NGC5533['blackhole']['t'], NGC5533['blackhole']['c'], NGC5533['blackhole']['k'])

# Bulge #############################
NGC5533['bulge'] = {
    'r' : np.asarray(NGC5533['raw_bulge']['xx']),
    'v' : np.asarray(NGC5533['raw_bulge']['yy']),
    'n'      : 2.7,     # Concentration parameter [unitless] (Noordermeer, Van Der Hulst, 2007)
    'q'      : 0.33,    # Intrinsic axis ratio [unitless] (Noordermeer, 2008)
    're_arcsec' : 9.9,     # Effective radius [arcsec] (Noordermeer, Van Der Hulst, 2007)
    'Lb'     : 3.6e10,  # Bulge luminosity [solar luminosity] (Calculated from absolute magnitude of the bulge)
    'ML'     : 2.8      # Mass-to-light ratio of bulge [unitless] (Noordermeer, 2008)
}
NGC5533['bulge']['t'], NGC5533['bulge']['c'], NGC5533['bulge']['k'] = inter.splrep(NGC5533['bulge']['r'], NGC5533['bulge']['v'])
NGC5533['bulge']['spline'] = inter.BSpline(NGC5533['bulge']['t'], NGC5533['bulge']['c'], NGC5533['bulge']['k'])
NGC5533['bulge']['re_rad'] = NGC5533['bulge']['re_arcsec'] * (np.pi / (3600*180)) #arcsec to rad
NGC5533['D_kpc'] = NGC5533['D_Mpc'] * 1000 #Mpc to kpc
NGC5533['bulge']['re_kpc'] = NGC5533['bulge']['re_rad'] * NGC5533['D_kpc'] #Effective radius kpc
NGC5533['sources'] = {
    'Mabs': ' [unitless] (Source #2, Table A4)',
    'n'   : ' [unitless] (Source #2, Table A4)',
    'q'   : ' [unitless] (Source #1, Table 1)',
    'i'   : ' [degrees] (Source #2, Table A2)',
    're'  : ' [kpc] (Source #2, Table A4)',
    'Lb'  : ' [Lsun] (Source #2, calculated from absolute magnitude)',
    'ML'  : ' [unitless] (Source #1, Table 1)',
    'D'   : ' [Mpc] (Source #2, Table 1)'
}

# Disk ##############################
NGC5533['disk'] = {
    'r' : np.asarray(NGC5533['raw_disk']['xx']),
    'v' : np.asarray(NGC5533['raw_disk']['yy'])
}
NGC5533['disk']['t'], NGC5533['disk']['c'], NGC5533['disk']['k'] = inter.splrep(NGC5533['disk']['r'], NGC5533['disk']['v'])
NGC5533['disk']['spline'] = inter.BSpline(NGC5533['disk']['t'], NGC5533['disk']['c'], NGC5533['disk']['k'])

# Halo ##############################
NGC5533['halo'] = {
    'r' : np.asarray(NGC5533['raw_halo']['xx']),
    'v' : np.asarray(NGC5533['raw_halo']['yy'])
}
NGC5533['halo']['t'], NGC5533['halo']['c'], NGC5533['halo']['k'] = inter.splrep(NGC5533['halo']['r'], NGC5533['halo']['v'])
NGC5533['halo']['spline'] = inter.BSpline(NGC5533['halo']['t'], NGC5533['halo']['c'], NGC5533['halo']['k'])

# Gas ###############################
NGC5533['gas'] = {
    'r' : np.asarray(NGC5533['raw_gas']['xx']),
    'v' : np.asarray(NGC5533['raw_gas']['yy'])
}
NGC5533['gas']['t'], NGC5533['gas']['c'], NGC5533['gas']['k'] = inter.splrep(NGC5533['gas']['r'], NGC5533['gas']['v'])
NGC5533['gas']['spline'] = inter.BSpline(NGC5533['gas']['t'], NGC5533['gas']['c'], NGC5533['gas']['k'])


###############################
########### NGC0891 ###########
###############################

NGC0891 = {
    # SPARC data files
    'raw_bulge'        : dp.getXYdata('data/NGC0891/0891_gBulge.dat'     ),   
    'raw_disk'         : dp.getXYdata('data/NGC0891/0891_gDisk.dat'      ),
    'raw_gas'          : dp.getXYdata('data/NGC0891/0891_gGas.dat'       ),

    # Get data
    'measured_data'    : dp.getXYdata_wYerr('data/NGC0891/0891_measured.dat'),
    
    # Some constants
    # Source: https://www.aanda.org/articles/aa/pdf/2011/07/aa16634-11.pdf 
    'Mabs'             : -21.66, # Absolute magnitude [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    'i'                : 89,     # Inclination angle [degrees] (Fraternali, Sancisi, and Kamphuis, 2011)
    'D_Mpc'            : 9.5,    # Distance [Mpc] (need this)
    'check'            : True
}
"""Data for galaxy NGC891. Parameters and data measurements in this dictionary are from [SPARC2016]_ unless noted otherwise.

:keys:   
    bulge: [dict] Further information pertaining to the bulge component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the bulge component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the bulge component.

        Lb: [float] Luminosity (solar luminosities) of the bulge.

        ML: [float] Mass-to-light ratio of the bulge.

        n: [float] Concentration parameter.

        q: [float] Intrinsic axis ratio.

        r: [array] Radial (kpc) values of the bulge component.

        re_arcsec: [float] Effective radius (arcsec).

        re_rad: [float] Effective radius (radians).

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the bulge component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the bulge component.

        v: [array] Velocity (km/s) values of the bulge component.

    check: [bool] A flag indicating that the necessary data for the "Literature Search Check" in `07_Bonus_Bulge_Rotation_Curve.ipynb <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/07_Bonus_Bulge_Rotation_Curve.ipynb>`_ is available.
    
    disk: [dict] Further information pertaining to the disk component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the disk component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the disk component.
        
        r: [array] Radial (kpc) values of the disk component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the disk component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the disk component.

        v: [array] Velocity (km/s) values of the disk component.

    D_Mpc: [float] The distance to the galaxy in Megaparsecs.
    
    galaxyname: [string] Name of the galaxy.

    gas: [dict] Further information pertaining to the gas component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the gas component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the gas component.
        
        r: [array] Radial (kpc) values of the gas component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the gas component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the gas component.

        v: [array] Velocity (km/s) values of the gas component.
    
    i: [int] The inclination angle, in degrees, of the galaxy.
    
    Mabs: [float] The absolute magnitude of the galaxy.

    massBH: [float] Mass of central black hole (solar masses).

    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).

    raw_bulge: [dict] Data representing the bulge's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    raw_disk: [dict] Data representing the disk's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.
    
    raw_gas: [dict] Data representing the gas' contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    rc: [float] Core radius (kpc).

    rho0: [float] Central mass density (solar mass/kpc^3).

    sources: [dict] Printing helper variables indicating units and sources for certain variables.
        D: [string] Units and source for distance to galaxy.
        
        i: [string] Units and source for inclination angle.

        Lb: [string] Units and source for bulge luminosity.
        
        Mabs: [string] Units and source for absolute magnitude.

        ML: [string] Units and source for bulge mass-light ratio.

        n: [string] Units and source for concentration parameter.

        re: [string] Units and source for effective radius.

        q: [string] Units and source for intrinsic axis ratio.
"""

# Parameters ########################
NGC0891['galaxyname'] = 'NGC 891'    # NGC catalog number of the galaxy
NGC0891['rho0'] = 1.4461e+08   # central mass density (in solar mass/kpc^3), calculated by fitting  
NGC0891['rc'] = 2.34           # core radius (in kpc), calculated by fitting 
NGC0891['massbh'] = 0          # central black hole is included in the bulge curve

#Organize measured data
NGC0891['m_radii']      = np.asarray(NGC0891['measured_data']['xx'])
NGC0891['m_velocities'] = np.asarray(NGC0891['measured_data']['yy'])
NGC0891['m_v_errors']   = np.asarray(NGC0891['measured_data']['ey'])

# Bulge #############################
NGC0891['bulge'] = {
    'r' : np.asarray(NGC0891['raw_bulge']['xx']),
    'v' : np.asarray(NGC0891['raw_bulge']['yy']),
    'n' : 2.99, # Concentration parameter [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    'q' : 0.68, # Intrinsic axis ratio [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    're_kpc' : 1.8,    # Effective radius [kpc] (Fraternali, Sancisi, and Kamphuis, 2011)
    'Lb'     : 2.2e10, # Bulge luminosity [solar luminosity] (Fraternali, Sancisi, and Kamphuis, 2011)
    'ML'     : 1.63,   # Mass-to-light ratio of bulge (Fraternali, Sancisi, and Kamphuis, 2011)
}
NGC0891['bulge']['t'], NGC0891['bulge']['c'], NGC0891['bulge']['k'] = inter.splrep(NGC0891['bulge']['r'], NGC0891['bulge']['v'])
NGC0891['bulge']['spline'] = inter.BSpline(NGC0891['bulge']['t'], NGC0891['bulge']['c'], NGC0891['bulge']['k'])
NGC0891['sources'] = {
    'Mabs': ' [unitless] (Source #3, calculated from Luminosity)',
    'n'   : ' [unitless] (Source #3, Table 4)',
    'q'   : ' [unitless] (Source #3, Table 4)',
    'i'   : ' [degrees] (Source #3, Table 1)',
    're'  : ' [kpc] (Source #3, Table 4)',
    'Lb'  : ' [Lsun] (Source #3, Table 4)',
    'ML'  : ' [unitless] (Source #3, Table 5, with dark matter halo)',
    'D'   : ' [Mpc] (Source #3, Table 1)'
}

# Disk ##############################
NGC0891['disk'] = {
    'r' : np.asarray(NGC0891['raw_disk']['xx']),
    'v' : np.asarray(NGC0891['raw_disk']['yy'])     
}
NGC0891['disk']['t'], NGC0891['disk']['c'], NGC0891['disk']['k'] = inter.splrep(NGC0891['disk']['r'], NGC0891['disk']['v'])
NGC0891['disk']['spline'] = inter.BSpline(NGC0891['disk']['t'], NGC0891['disk']['c'], NGC0891['disk']['k'])

# Gas ###############################
NGC0891['gas'] = {
    'r' : np.asarray(NGC0891['raw_gas']['xx']),
    'v' : np.asarray(NGC0891['raw_gas']['yy'])
}
NGC0891['gas']['t'], NGC0891['gas']['c'], NGC0891['gas']['k'] = inter.splrep(NGC0891['gas']['r'], NGC0891['gas']['v'])
NGC0891['gas']['spline'] = inter.BSpline(NGC0891['gas']['t'], NGC0891['gas']['c'], NGC0891['gas']['k'])

NGC891 = NGC0891    # Considering when someone doesn't type the 0
"""Alias for :func:`NGC0891 <load_galaxies.NGC0891>`.
"""

###############################
########### NGC7814 ###########
###############################

NGC7814 = {
    # SPARC data files
    'raw_bulge'        : dp.getXYdata('data/NGC7814/7814_gBulge.dat'     ),
    'raw_disk'         : dp.getXYdata('data/NGC7814/7814_gDisk.dat'      ),
    'raw_gas'          : dp.getXYdata('data/NGC7814/7814_gGas.dat'       ),

    # Get data
    'measured_data'    : dp.getXYdata_wYerr('data/NGC7814/7814_measured.dat'),
    
    # Some constants
    # Source: https://www.aanda.org/articles/aa/pdf/2011/07/aa16634-11.pdf 
    'i'                : 90,   # Inclination angle [degrees] (Fraternali, Sancisi, and Kamphuis, 2011)
    'D_Mpc'            : 14.6,  # Distance [Mpc]
    'Mabs'             : -22.38,# Absolute magnitude [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    'check'            : True
}
"""Data for galaxy NGC7814. Parameters and data measurements in this dictionary are from [SPARC2016]_ unless noted otherwise.

:keys:
    bulge: [dict] Further information pertaining to the bulge component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the bulge component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the bulge component.

        Lb: [float] Luminosity (solar luminosities) of the bulge.

        ML: [float] Mass-to-light ratio of the bulge.

        n: [float] Concentration parameter.

        q: [float] Intrinsic axis ratio.

        r: [array] Radial (kpc) values of the bulge component.

        re_arcsec: [float] Effective radius (arcsec).

        re_rad: [float] Effective radius (radians).

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the bulge component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the bulge component.

        v: [array] Velocity (km/s) values of the bulge component.

    check: [bool] A flag indicating that the necessary data for the "Literature Search Check" in `07_Bonus_Bulge_Rotation_Curve.ipynb <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/07_Bonus_Bulge_Rotation_Curve.ipynb>`_ is available.
    
    disk: [dict] Further information pertaining to the disk component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the disk component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the disk component.
        
        r: [array] Radial (kpc) values of the disk component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the disk component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the disk component.

        v: [array] Velocity (km/s) values of the disk component.

    D_Mpc: [float] The distance to the galaxy in Megaparsecs.
    
    galaxyname: [string] Name of the galaxy.

    gas: [dict] Further information pertaining to the gas component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the gas component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the gas component.
        
        r: [array] Radial (kpc) values of the gas component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the gas component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the gas component.

        v: [array] Velocity (km/s) values of the gas component.
    
    i: [int] The inclination angle, in degrees, of the galaxy.
    
    Mabs: [float] The absolute magnitude of the galaxy.

    massBH: [float] Mass of central black hole (solar masses).

    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).

    raw_bulge: [dict] Data representing the bulge's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    raw_disk: [dict] Data representing the disk's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.
    
    raw_gas: [dict] Data representing the gas' contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    rc: [float] Core radius (kpc).

    rho0: [float] Central mass density (solar mass/kpc^3).

    sources: [dict] Printing helper variables indicating units and sources for certain variables.
        D: [string] Units and source for distance to galaxy.
        
        i: [string] Units and source for inclination angle.

        Lb: [string] Units and source for bulge luminosity.
        
        Mabs: [string] Units and source for absolute magnitude.

        ML: [string] Units and source for bulge mass-light ratio.

        n: [string] Units and source for concentration parameter.

        re: [string] Units and source for effective radius.

        q: [string] Units and source for intrinsic axis ratio.
"""

# Parameters ########################
NGC7814['galaxyname'] = 'NGC 7814'   # NGC catalog number of the galaxy
NGC7814['rho0'] = 1.524e8       # central mass density (in solar mass/kpc^3), Source: Fraternali et al. (2011) 
NGC7814['rc'] = 2.1            # core radius (in kpc), Source: Fraternali et al. (2011) 
NGC7814['massbh'] = 0          # central black hole is included in the bulge curve

#Organize measured data
NGC7814['m_radii']      = np.asarray(NGC7814['measured_data']['xx'])
NGC7814['m_velocities'] = np.asarray(NGC7814['measured_data']['yy'])
NGC7814['m_v_errors']   = np.asarray(NGC7814['measured_data']['ey'])

# Bulge #############################
NGC7814['bulge'] = {
    'r'      : np.asarray(NGC7814['raw_bulge']['xx']),
    'v'      : np.asarray(NGC7814['raw_bulge']['yy']),
    'n'      : 4.0,      # Concentration parameter [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    'q'      : 0.61,     # Intrinsic axis ratio [unitless] (Fraternali, Sancisi, and Kamphuis, 2011)
    're_kpc' : 2.16,     # Effective radius [kpc] (Fraternali, Sancisi, and Kamphuis, 2011)
    'Lb'     : 7e10,     # Bulge luminosity [solar luminosity] (Fraternali, Sancisi, and Kamphuis, 2011)
    'ML'     : 0.71      # Mass-to-light ratio of bulge (Fraternali, Sancisi, and Kamphuis, 2011)
}
NGC7814['bulge']['t'], NGC7814['bulge']['c'], NGC7814['bulge']['k'] = inter.splrep(NGC7814['bulge']['r'], NGC7814['bulge']['v'])
NGC7814['bulge']['spline'] = inter.BSpline(NGC7814['bulge']['t'], NGC7814['bulge']['c'], NGC7814['bulge']['k'])
NGC7814['sources'] = NGC0891['sources'] #The entries were identical.

# Disk ##############################
NGC7814['disk'] = {
    'r' : np.asarray(NGC7814['raw_disk']['xx']),
    'v' : np.asarray(NGC7814['raw_disk']['yy'])
}
NGC7814['disk']['t'], NGC7814['disk']['c'], NGC7814['disk']['k'] = inter.splrep(NGC7814['disk']['r'], NGC7814['disk']['v'])
NGC7814['disk']['spline'] = inter.BSpline(NGC7814['disk']['t'], NGC7814['disk']['c'], NGC7814['disk']['k'])

# Gas ###############################
NGC7814['gas'] = {
    'r' : np.asarray(NGC7814['raw_gas']['xx']),
    'v' : np.asarray(NGC7814['raw_gas']['yy'])
}
NGC7814['gas']['t'], NGC7814['gas']['c'], NGC7814['gas']['k'] = inter.splrep(NGC7814['gas']['r'], NGC7814['gas']['v'])
NGC7814['gas']['spline'] = inter.BSpline(NGC7814['gas']['t'], NGC7814['gas']['c'], NGC7814['gas']['k'])


###############################
########### NGC5005 ###########
###############################

NGC5005 = {
    
    'raw_bulge'        : dp.getXYdata('data/NGC5005/ngc5005_bulge.txt'     ),
    'raw_disk'         : dp.getXYdata('data/NGC5005/ngc5005_disk.txt'      ),
    'raw_halo'         : dp.getXYdata('data/NGC5005/ngc5005_halo.txt'      ),
    'raw_gas'          : dp.getXYdata('data/NGC5005/ngc5005_gas.txt'       ),

    # Get data
    'measured_data'    : dp.getXYdata_wXYerr('data/NGC5005/ngc5005_data.txt')
}
"""Data for galaxy NGC5005. Parameters and data measurements in this dictionary are from [SPARC2016]_ unless noted otherwise.

:keys:
    bulge: [dict] Further information pertaining to the bulge component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the bulge component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the bulge component.

        Lb: [float] Luminosity (solar luminosities) of the bulge.

        ML: [float] Mass-to-light ratio of the bulge.

        n: [float] Concentration parameter.

        q: [float] Intrinsic axis ratio.

        r: [array] Radial (kpc) values of the bulge component.

        re_arcsec: [float] Effective radius (arcsec).

        re_rad: [float] Effective radius (radians).

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the bulge component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the bulge component.

        v: [array] Velocity (km/s) values of the bulge component.

    check: [bool] A flag indicating that the necessary data for the "Literature Search Check" in `07_Bonus_Bulge_Rotation_Curve.ipynb <https://github.com/villano-lab/galactic-spin-W1/blob/master/binder/07_Bonus_Bulge_Rotation_Curve.ipynb>`_ is available.
    
    disk: [dict] Further information pertaining to the disk component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the disk component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the disk component.
        
        r: [array] Radial (kpc) values of the disk component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the disk component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the disk component.

        v: [array] Velocity (km/s) values of the disk component.

    D_Mpc: [float] The distance to the galaxy in Megaparsecs.

    halo: [dict] Further information pertaining to the halo component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the halo component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the halo component.
        
        r: [array] Radial (kpc) values of the halo component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the halo component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the halo component.

        v: [array] Velocity (km/s) values of the halo component.
    
    galaxyname: [string] Name of the galaxy.

    gas: [dict] Further information pertaining to the gas component.
        c: [array] B-spline coefficients returned by using scipy.interpolate.splrep on the gas component.

        k: [int] Degree of spline returned by using scipy.interpolate.splrep on the gas component.
        
        r: [array] Radial (kpc) values of the gas component.

        spline: [scipy.interpolate._bsplines.BSpline] Spline of the gas component.
        
        t: [array] Vector of knots returned by using scipy.interpolate.splrep on the gas component.

        v: [array] Velocity (km/s) values of the gas component.
    
    i: [int] The inclination angle, in degrees, of the galaxy.
    
    Mabs: [float] The absolute magnitude of the galaxy.

    massBH: [float] Mass of central black hole (solar masses).

    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).

    raw_bulge: [dict] Data representing the bulge's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    raw_disk: [dict] Data representing the disk's contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.
    
    raw_gas: [dict] Data representing the gas' contribution to the theoretical rotation curve.  
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

    rc: [float] Core radius (kpc).

    rho0: [float] Central mass density (solar mass/kpc^3).
"""

#Organize measured data
NGC5005['m_radii']      = np.asarray(NGC5005['measured_data']['xx'])
NGC5005['m_velocities'] = np.asarray(NGC5005['measured_data']['yy'])
NGC5005['m_r_errors']   = np.asarray(NGC5005['measured_data']['ex'])
NGC5005['m_v_errors']   = np.asarray(NGC5005['measured_data']['ey'])

# Bulge #############################
NGC5005['bulge'] = {
    'r' : np.asarray(NGC5005['raw_bulge']['xx']),
    'v' : np.asarray(NGC5005['raw_bulge']['yy'])
}
NGC5005['bulge']['t'], NGC5005['bulge']['c'], NGC5005['bulge']['k'] = inter.splrep(NGC5005['bulge']['r'], NGC5005['bulge']['v'])
NGC5005['bulge']['spline'] = inter.BSpline(NGC5005['bulge']['t'], NGC5005['bulge']['c'], NGC5005['bulge']['k'])

# Disk ##############################
NGC5005['disk'] = {
    'r' : np.asarray(NGC5005['raw_disk']['xx']),
    'v' : np.asarray(NGC5005['raw_disk']['yy'])
}
NGC5005['disk']['t'], NGC5005['disk']['c'], NGC5005['disk']['k'] = inter.splrep(NGC5005['disk']['r'], NGC5005['disk']['v'])
NGC5005['disk']['spline'] = inter.BSpline(NGC5005['disk']['t'], NGC5005['disk']['c'], NGC5005['disk']['k'])

# Halo ##############################
NGC5005['halo'] = {
    'r' : np.asarray(NGC5005['raw_halo']['xx']),
    'v' : np.asarray(NGC5005['raw_halo']['yy'])
}
NGC5005['halo']['t'], NGC5005['halo']['c'], NGC5005['halo']['k'] = inter.splrep(NGC5005['halo']['r'], NGC5005['halo']['v'])
NGC5005['halo']['spline'] = inter.BSpline(NGC5005['halo']['t'], NGC5005['halo']['c'], NGC5005['halo']['k'])

# Gas ###############################
NGC5005['gas'] = {
    'r' : np.asarray(NGC5005['raw_gas']['xx']),
    'v' : np.asarray(NGC5005['raw_gas']['yy'])
}
NGC5005['gas']['t'], NGC5005['gas']['c'], NGC5005['gas']['k'] = inter.splrep(NGC5005['gas']['r'], NGC5005['gas']['v'])
NGC5005['gas']['spline'] = inter.BSpline(NGC5005['gas']['t'], NGC5005['gas']['c'], NGC5005['gas']['k'])

# Parameters ########################
NGC5005['galaxyname'] = 'NGC 5005'   # NGC catalog number of the galaxy
NGC5005['rho0'] = 1.693e+07          # central mass density (in solar mass/kpc^3), obtained by fitting   
NGC5005['rc'] = 9.917                # core radius (in kpc), obtained by fitting 
NGC5005['massbh'] = 0                # central black hole is included in the bulge curve


###############################
####### Other Galaxies ########
###############################

# NGC 3198
NGC3198 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/NGC3198.txt')}
"""Data for galaxy NGC3198. Parameters and data measurements in this dictionary are from [Karukes2015]_ unless noted otherwise.

:keys:
    galaxyname: [string] Name of the galaxy.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).
    
    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.
    
"""
NGC3198['m_radii']      = np.asarray(NGC3198['measured_data']['xx'])
NGC3198['m_velocities'] = np.asarray(NGC3198['measured_data']['yy'])
NGC3198['m_v_errors']   = np.asarray(NGC3198['measured_data']['ey'])
NGC3198['galaxyname'] = 'NGC 3198' 


# UGC 477
UGC0477 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC477.txt')}
"""Data for galaxy UGC477. Parameters and data measurements in this dictionary are from [deNaray2006]_ and [deNaray2008]_ unless noted otherwise.

:keys:
    galaxyname: [string] Name of the galaxy.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).
    
    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.
    
"""
UGC0477['m_radii']      = np.asarray(UGC0477['measured_data']['xx'])
UGC0477['m_velocities'] = np.asarray(UGC0477['measured_data']['yy'])
UGC0477['m_v_errors']   = np.asarray(UGC0477['measured_data']['ey'])
UGC0477['galaxyname'] = 'UGC 477'
UGC477 = UGC0477
"""Alias for :func:`UGC0477 <load_galaxies.UGC0477>`.
"""

# UGC 1281
UGC1281 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC1281.txt')}
"""Data for galaxy UGC1281. Parameters and data measurements in this dictionary are from [deNaray2006]_ and [deNaray2008]_ unless noted otherwise.

:keys:
    galaxyname: [string] Name of the galaxy.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).
    
    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.
    
"""
UGC1281['m_radii']      = np.asarray(UGC1281['measured_data']['xx'])
UGC1281['m_velocities'] = np.asarray(UGC1281['measured_data']['yy'])
UGC1281['m_v_errors']   = np.asarray(UGC1281['measured_data']['ey'])
UGC1281['galaxyname'] = 'UGC 1281' 

# UGC 1437
UGC1437 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC1437.txt')}
"""Data for galaxy UGC1437. Parameters and data measurements in this dictionary are from [Epinat2008]_ unless noted otherwise.

:keys:
    galaxyname: [string] Name of the galaxy.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).
    
    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.
    
"""
UGC1437['m_radii']      = np.asarray(UGC1437['measured_data']['xx'])
UGC1437['m_velocities'] = np.asarray(UGC1437['measured_data']['yy'])
UGC1437['m_v_errors']   = np.asarray(UGC1437['measured_data']['ey'])
UGC1437['galaxyname'] = 'UGC 1437' 

# UGC 2953
UGC2953 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC2953.txt')}
"""Data for galaxy UGC2953. Parameters and data measurements in this dictionary are from [SPARC2016]_ unless noted otherwise.

:keys:
    galaxyname: [string] Name of the galaxy.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).
    
    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.
    
"""
UGC2953['m_radii']      = np.asarray(UGC2953['measured_data']['xx'])
UGC2953['m_velocities'] = np.asarray(UGC2953['measured_data']['yy'])
UGC2953['m_v_errors']   = np.asarray(UGC2953['measured_data']['ey'])
UGC2953['galaxyname'] = 'UGC 2953'

# UGC 4325
UGC4325 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC4325.txt')}
"""Data for galaxy UGC4325. Parameters and data measurements in this dictionary are from [deNaray2006]_ and [deNaray2008]_ unless noted otherwise.

:keys:
    galaxyname: [string] Name of the galaxy.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).
    
    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.
    
"""
UGC4325['m_radii']      = np.asarray(UGC4325['measured_data']['xx'])
UGC4325['m_velocities'] = np.asarray(UGC4325['measured_data']['yy'])
UGC4325['m_v_errors']   = np.asarray(UGC4325['measured_data']['ey'])
UGC4325['galaxyname'] = 'UGC 4325'

# UGC 5253
UGC5253 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC5253.txt')}
"""Data for galaxy UGC5253. Parameters and data measurements in this dictionary are from [SPARC2016]_ unless noted otherwise.

:keys:
    galaxyname: [string] Name of the galaxy.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).
    
    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.
    
"""
UGC5253['m_radii']      = np.asarray(UGC5253['measured_data']['xx'])
UGC5253['m_velocities'] = np.asarray(UGC5253['measured_data']['yy'])
UGC5253['m_v_errors']   = np.asarray(UGC5253['measured_data']['ey'])
UGC5253['galaxyname'] = 'UGC 5253'   

# UGC 6787
UGC6787 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC6787.txt')}
"""Data for galaxy UGC6787. Parameters and data measurements in this dictionary are from [SPARC2016]_ unless noted otherwise.

:keys:
    galaxyname: [string] Name of the galaxy.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).
    
    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.
    
"""
UGC6787['m_radii']      = np.asarray(UGC6787['measured_data']['xx'])
UGC6787['m_velocities'] = np.asarray(UGC6787['measured_data']['yy'])
UGC6787['m_v_errors']   = np.asarray(UGC6787['measured_data']['ey'])
UGC6787['galaxyname'] = 'UGC 6787'   

# UGC 10075
UGC10075 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC10075.txt')}
"""Data for galaxy UGC10075. Parameters and data measurements in this dictionary are from [Epinat2008]_.

:keys:
    galaxyname: [string] Name of the galaxy.

    m_radii: [array] Radii of measured data (kpc).

    m_v_errors: [array] Errors of measured velocities (km/s).

    m_velocities: [array] Velocities of measured data (km/s).
    
    measured_data: [dict] Data representing the bottom of the error band on the total curve.
        xx: [list] Float values representing the radial (kpc) data.  

        yy: [list] Float values representing the velocity (km/s) data.

        ey: [list] Float values representing the error on the radial (kpc) data.
    
"""
UGC10075['m_radii']      = np.asarray(UGC10075['measured_data']['xx'])
UGC10075['m_velocities'] = np.asarray(UGC10075['measured_data']['yy'])
UGC10075['m_v_errors']   = np.asarray(UGC10075['measured_data']['ey'])
UGC10075['galaxyname'] = 'UGC 10075'