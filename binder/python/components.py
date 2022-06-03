###############
### Imports ###
###############

import numpy as np
import dataPython as dp
import scipy.integrate as si
from scipy.interpolate import InterpolatedUnivariateSpline      # Spline function
import lmfit as lm                                              # Fitting

ngc7814list = ['ngc7814', 'ngc 7814']
ngc5533list = ['ngc5533', 'ngc 5533']

##############################
### Import data/text files ###
##############################

# Note: there's no need to import the radius for each component as everything has the same r array 
# (the r array of the raw data)

# Datapoints:
def data(galaxy):
    if galaxy.lower() in ngc7814list:
        return dp.getXYdata_wYerr('data/NGC7814/7814_measured.dat')
    elif galaxy.lower() in ngc5533list:
        return dataNGC5533.data
def r_dat(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(data(galaxy)['xx'])
    elif galaxy.lower() in ngc5533list:
        return dataNGC5533.r_dat
def v_dat(galaxy):
    if galaxy.lower() in ngc7814list:
        np.asarray(data(galaxy)['yy'])
    elif galaxy.lower() in ngc5533list:
def v_err0(galaxy):
    if galaxy.lower() in ngc5533list:
        return dataNGC5533.v_err0
    else:
        return None
def v_err1(galaxy):
    if galaxy.lower() in ngc7814list:
        return 
        