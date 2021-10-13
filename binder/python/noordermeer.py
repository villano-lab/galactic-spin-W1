# Traced curves from Noordermeer's paper: Rotation curves of flattened Sersic bulges, Figure 4.

################################
########### Imports ############
################################
import sys
sys.path.append('../../python/')

import dataPython as dp
import numpy as np
import scipy.interpolate as inter
import matplotlib.pyplot as plt

################################
######### Data files ###########
################################
data_total = dp.getXYdata('../data/final/nord-120kpc-total.txt')
data_bh = dp.getXYdata('../data/final/nord-120kpc-blackhole.txt')
data_bulge = dp.getXYdata('../data/final/nord-120kpc-bulge.txt')
data_disk = dp.getXYdata('../data/final/nord-120kpc-disk.txt')
data_halo = dp.getXYdata('../data/final/nord-120kpc-halo.txt')
data_gas = dp.getXYdata('../data/final/nord-120kpc-gas.txt')
data_greyb_bottom = dp.getXYdata('../data/final/nord-120kpc-bottomband.txt')
data_greyb_top = dp.getXYdata('../data/final/nord-120kpc-topband.txt')

rval = np.linspace(0.1,120,500)
rb = np.linspace(0.1,120,500)

################################
##### Measured data points #####
################################
data = dp.getXYdata_wXYerr('../data/100kpc_data.txt')
r_dat = np.asarray(data['xx'])
v_dat = np.asarray(data['yy'])
v_err0 = np.asarray(data['ex'])
v_err1 = np.asarray(data['ey'])

################################
####### Uncertainty band #######
################################
# convert to numpy arrays
r_bottomband = np.asarray(data_greyb_bottom['xx'])
v_bottomband = np.asarray(data_greyb_bottom['yy'])
r_topband = np.asarray(data_greyb_top['xx'])
v_topband = np.asarray(data_greyb_top['yy'])

band = (v_topband - v_bottomband)/2
# for weightdata, lengths of v_err1 and band must equal, make band array same length as v_err1
band = band[0::28]
band = band[1:]

# smoothing - new, `spline` would not run on my computer
tb, cb, kb = inter.splrep(r_bottomband,v_bottomband)
tt, ct, kt = inter.splrep(r_topband,v_topband)

greyb_bottom = inter.BSpline(tb, cb, kb)
greyb_top    = inter.BSpline(tt, ct, kt)

################################
######### Total curve ##########
################################
r_total = np.asarray(data_total['xx'])
v_total = np.asarray(data_total['yy'])

tt, ct, kt = inter.splrep(r_total,v_total)
noord_total = inter.BSpline(tt,ct,kt)

################################
######### Black Hole ###########
################################
r_bh = np.asarray(data_bh['xx'])
v_bh = np.asarray(data_bh['yy'])

tbh, cbh, kbh = inter.splrep(r_bh,v_bh)
noord_blackhole = inter.BSpline(tbh,cbh,kbh)

################################
############ Bulge #############
################################
r_bulge = np.asarray(data_bulge['xx'])
v_bulge = np.asarray(data_bulge['yy'])

tb, cb, kb = inter.splrep(r_bulge,v_bulge)
noord_bulge = inter.BSpline(tb,cb,kb)

################################
############ Disk ##############
################################
r_disk = np.asarray(data_disk['xx'])
v_disk = np.asarray(data_disk['yy'])

td, cd, kd = inter.splrep(r_disk,v_disk)
noord_disk = inter.BSpline(td,cd,kd)

################################
############ Halo ##############
################################
r_halo = np.asarray(data_halo['xx'])
v_halo = np.asarray(data_halo['yy'])

th, ch, kh = inter.splrep(r_halo,v_halo)
noord_halo = inter.BSpline(th,ch,kh)

################################
############# Gas ##############
################################
r_gas = np.asarray(data_gas['xx'])
v_gas = np.asarray(data_gas['yy'])

tg,cg,kg = inter.splrep(r_gas,v_gas)
noord_gas = inter.BSpline(tg,cg,kg)
