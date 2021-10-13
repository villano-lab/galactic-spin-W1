#Data to fit to for each galaxy to be used in workshop

###############################
########## Imports ############
###############################
import sys
sys.path.append('../python/')

import dataPython         as dp
import numpy              as np
import scipy.interpolate  as inter
import matplotlib.pyplot  as plt


###############################
########### NGC5533 ###########
###############################

NGC5533 = {

    # Load data from files for Noordermeer's band and fitted curves
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
    'measured_data'    : dp.getXYdata_wXYerr('data/NGC5533/100kpc_data.txt')
        
}

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
NGC5533['n_v_bandwidth'] = NGC5533['n_v_topband'] - NGC5533['n_v_btmband']
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
    'r' : np.asarray(NGC5533['raw_blackhole']['xx']),
    'v' : np.asarray(NGC5533['raw_blackhole']['yy'])
}
NGC5533['blackhole']['t'], NGC5533['blackhole']['c'], NGC5533['blackhole']['k'] = inter.splrep(NGC5533['blackhole']['r'], NGC5533['blackhole']['v'])
NGC5533['blackhole']['spline'] = inter.BSpline(NGC5533['blackhole']['t'], NGC5533['blackhole']['c'], NGC5533['blackhole']['k'])

# Bulge #############################
NGC5533['bulge'] = {
    'r' : np.asarray(NGC5533['raw_bulge']['xx']),
    'v' : np.asarray(NGC5533['raw_bulge']['yy'])
}
NGC5533['bulge']['t'], NGC5533['bulge']['c'], NGC5533['bulge']['k'] = inter.splrep(NGC5533['bulge']['r'], NGC5533['bulge']['v'])
NGC5533['bulge']['spline'] = inter.BSpline(NGC5533['bulge']['t'], NGC5533['bulge']['c'], NGC5533['bulge']['k'])

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
    
    'raw_bulge'        : dp.getXYdata('data/NGC0891/891_dtBulge.dat'     ),
    'raw_disk'         : dp.getXYdata('data/NGC0891/891_dtDisk.dat'      ),
    'raw_gas'          : dp.getXYdata('data/NGC0891/891_dtGas.dat'       ),

    # Get data
    'measured_data'    : dp.getXYdata_wXYerr('data/NGC0891/891_data')
}

# Parameters ########################
NGC0891['galaxyname'] = 'NGC 891'    # NGC catalog number of the galaxy
NGC0891['rho0'] = 3.31e7       # central mass density (in solar mass/kpc^3), Source: Richards et al. (2015)    
NGC0891['rc'] = 1.9            # core radius (in kpc), Source: Richards et al. (2015)
NGC0891['massbh'] = 0          # central black hole is included in the bulge curve

#Organize measured data
NGC0891['m_radii']      = np.asarray(NGC0891['measured_data']['xx'])
NGC0891['m_velocities'] = np.asarray(NGC0891['measured_data']['yy'])
NGC0891['m_r_errors']   = np.asarray(NGC0891['measured_data']['ex'])
NGC0891['m_v_errors']   = np.asarray(NGC0891['measured_data']['ey'])

# Bulge #############################
NGC0891['bulge'] = {
    'r' : np.asarray(NGC0891['raw_bulge']['xx']),
    'v' : np.asarray(NGC0891['raw_bulge']['yy'])
}
NGC0891['bulge']['t'], NGC0891['bulge']['c'], NGC0891['bulge']['k'] = inter.splrep(NGC0891['bulge']['r'], NGC0891['bulge']['v'])
NGC0891['bulge']['spline'] = inter.BSpline(NGC0891['bulge']['t'], NGC0891['bulge']['c'], NGC0891['bulge']['k'])

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

NGC891 = NGC0891    # Considering when someone forgets to type 0

###############################
########### NGC7814 ###########
###############################

NGC7814 = {
    
    'raw_bulge'        : dp.getXYdata('data/NGC7814/7814reallybulge.dat'     ),
    'raw_disk'         : dp.getXYdata('data/NGC7814/7814reallydisk.dat'      ),
    'raw_gas'          : dp.getXYdata('data/NGC7814/7814reallygas.dat'       ),

    # Get data
    'measured_data'    : dp.getXYdata_wXYerr('data/NGC7814/ngc7814data')
}

# Parameters ########################
NGC7814['galaxyname'] = 'NGC 7814'   # NGC catalog number of the galaxy
NGC7814['rho0'] = 1.52e8       # central mass density (in solar mass/kpc^3), Source: Richards et al. (2015)
NGC7814['rc'] = 2.1            # core radius (in kpc), Source: Richards et al. (2015)
NGC7814['massbh'] = 0          # central black hole is included in the bulge curve

#Organize measured data
NGC7814['m_radii']      = np.asarray(NGC7814['measured_data']['xx'])
NGC7814['m_velocities'] = np.asarray(NGC7814['measured_data']['yy'])
NGC7814['m_r_errors']   = np.asarray(NGC7814['measured_data']['ex'])
NGC7814['m_v_errors']   = np.asarray(NGC7814['measured_data']['ey'])

# Bulge #############################
NGC7814['bulge'] = {
    'r' : np.asarray(NGC7814['raw_bulge']['xx']),
    'v' : np.asarray(NGC7814['raw_bulge']['yy'])
}
NGC7814['bulge']['t'], NGC7814['bulge']['c'], NGC7814['bulge']['k'] = inter.splrep(NGC7814['bulge']['r'], NGC7814['bulge']['v'])
NGC7814['bulge']['spline'] = inter.BSpline(NGC7814['bulge']['t'], NGC7814['bulge']['c'], NGC7814['bulge']['k'])

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
NGC5005['rho0'] = 1e8       # central mass density (in solar mass/kpc^3), guess!    
NGC5005['rc'] = 2.5         # core radius (in kpc), Source: Richards et al. (2015)
NGC5005['massbh'] = 0       # central black hole is included in the bulge curve


###############################
####### Other Galaxies ########
###############################

# NGC 3198
NGC3198 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/NGC3198.txt')}
NGC3198['m_radii']      = np.asarray(NGC3198['measured_data']['xx'])
NGC3198['m_velocities'] = np.asarray(NGC3198['measured_data']['yy'])
NGC3198['m_v_errors']   = np.asarray(NGC3198['measured_data']['ey'])
NGC3198['galaxyname'] = 'NGC 3198' 

# UGC 89
#UGC89 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC89.txt')}
#UGC89['m_radii']      = np.asarray(UGC89['measured_data']['xx'])
#UGC89['m_velocities'] = np.asarray(UGC89['measured_data']['yy'])
#UGC89['m_v_errors']   = np.asarray(UGC89['measured_data']['ey'])
#UGC89['galaxyname'] = 'UGC 89' 

# UGC 477
UGC477 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC477.txt')}
UGC477['m_radii']      = np.asarray(UGC477['measured_data']['xx'])
UGC477['m_velocities'] = np.asarray(UGC477['measured_data']['yy'])
UGC477['m_v_errors']   = np.asarray(UGC477['measured_data']['ey'])
UGC477['galaxyname'] = 'UGC 477'

# UGC 1281
UGC1281 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC1281.txt')}
UGC1281['m_radii']      = np.asarray(UGC1281['measured_data']['xx'])
UGC1281['m_velocities'] = np.asarray(UGC1281['measured_data']['yy'])
UGC1281['m_v_errors']   = np.asarray(UGC1281['measured_data']['ey'])
UGC1281['galaxyname'] = 'UGC 1281' 

# UGC 1437
UGC1437 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC1437.txt')}
UGC1437['m_radii']      = np.asarray(UGC1437['measured_data']['xx'])
UGC1437['m_velocities'] = np.asarray(UGC1437['measured_data']['yy'])
UGC1437['m_v_errors']   = np.asarray(UGC1437['measured_data']['ey'])
UGC1437['galaxyname'] = 'UGC 1437' 

# UGC 2953
UGC2953 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC2953.txt')}
UGC2953['m_radii']      = np.asarray(UGC2953['measured_data']['xx'])
UGC2953['m_velocities'] = np.asarray(UGC2953['measured_data']['yy'])
UGC2953['m_v_errors']   = np.asarray(UGC2953['measured_data']['ey'])
UGC2953['galaxyname'] = 'UGC 2953'

# UGC 4325
UGC4325 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC4325.txt')}
UGC4325['m_radii']      = np.asarray(UGC4325['measured_data']['xx'])
UGC4325['m_velocities'] = np.asarray(UGC4325['measured_data']['yy'])
UGC4325['m_v_errors']   = np.asarray(UGC4325['measured_data']['ey'])
UGC4325['galaxyname'] = 'UGC 4325'

# UGC 5253
UGC5253 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC5253.txt')}
UGC5253['m_radii']      = np.asarray(UGC5253['measured_data']['xx'])
UGC5253['m_velocities'] = np.asarray(UGC5253['measured_data']['yy'])
UGC5253['m_v_errors']   = np.asarray(UGC5253['measured_data']['ey'])
UGC5253['galaxyname'] = 'UGC 5253'   

# UGC 6787
UGC6787 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC6787.txt')}
UGC6787['m_radii']      = np.asarray(UGC6787['measured_data']['xx'])
UGC6787['m_velocities'] = np.asarray(UGC6787['measured_data']['yy'])
UGC6787['m_v_errors']   = np.asarray(UGC6787['measured_data']['ey'])
UGC6787['galaxyname'] = 'UGC 6787'   

# UGC 10075
UGC10075 = {'measured_data' : dp.getXYdata_wYerr('data/othergalaxies/UGC10075.txt')}
UGC10075['m_radii']      = np.asarray(UGC10075['measured_data']['xx'])
UGC10075['m_velocities'] = np.asarray(UGC10075['measured_data']['yy'])
UGC10075['m_v_errors']   = np.asarray(UGC10075['measured_data']['ey'])
UGC10075['galaxyname'] = 'UGC 10075'