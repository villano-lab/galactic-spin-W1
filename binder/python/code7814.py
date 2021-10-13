#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm
import sys
sys.path.append('../../python/')
import dataPython as dp
import NGC5533_functions as nf
from sympy import *






#kpc limits, visual guess based on the galaxy in the image chosen:
minkpc=0
maxkpc=25


w=3970 #width of the square image
h=w
#center of galaxy:
c_x=w/2
c_y=h/2.2

#galaxy parameters
rcut= 2.1 #cutoff radius, kpc for ngc7814...
#rho00=154.4/1000 /(1e3**3) #Msun/kpc^3 for ngc7814 NOT CURRENTLY BEING USED. ALSO WHY IS IT SO SMALL???

GG = 4.30091e-6    #gravitational constant (kpc/solar mass*(km/s)^2)

#visual scaling
scale=23e6 #the number of theoretical black holes each graphed dot actually presents, somewhat arbitrary but 
#should be a number that does't require the plotting too many or too few dots representing bh's
#units: scale = [#number of actual black holes / plotted dot]
kpctopixels=50 #visual scaling, varies depending on size of galaxy image (and actual size of galaxy)
r1=minkpc*kpctopixels
r2=maxkpc*kpctopixels

#for number of black holes slider
Max=100 #max number of blackholes
best_A=.5*Max #default # of bh's for slider
stepN=.1*Max #step of # of bh's slider
Min=1# min # of black holes

#for blackhole mass slider:
minmass=.1 #solar masses, arbitrary
maxmass=3.7 #solar masses, just smaller then the smallest black hole ever discovered according to
#https://www.scientificamerican.com/gallery/the-smallest-known-black-hole/
start=.5*maxmass #default mass value for slider


                                  
#calculate the mass distribution on the galaxy using parametrs of similar galaxies (???)
#r = symbols('r')
#X=r/rcut
#Z = rho00/(X*(1+X)**2) #NFW (dark halo) density profile. I'm assuming Rs is cutoff radius. is this correct?
#RHO=4*np.pi*r**2*Z #
#RRR=integrate(RHO)
#func = lambdify(r, RRR,'numpy') #returns a numpy-ready function
data = dp.getXYdata_wXYerr('../testing/7814/ngc7814data')
r = np.asarray(data['xx'])
viso=(1-((rcut/r)*np.arctan(r/rcut))) #this is just the unitless part of the h_viso function in our library
rr=np.linspace(minkpc,maxkpc,int(Max)) #kpc
rrr = np.random.uniform(r1,r2,int(Max))

angle=np.random.uniform(0,2*np.pi,int(Max)) #angle 0 to 360 degrees for fulle circle (donut) for each bracket






# In[2]:


#**********************importing text files******************************
#there's no need to import the radius for each component as everything has the same r array (the r array of the raw data)

v_dat = np.asarray(data['yy'])
v_err1 = np.asarray(data['ey'])
#disk:
disk = dp.getXYZdata('../testing/7814/7814reallydisk.dat')
d = np.asarray(disk['zz'])
#bulge:
bulge = dp.getXYZdata('../testing/7814/7814reallybulge.dat')
b = np.asarray(bulge['zz'])
#gas:
gas = dp.getXYZdata('../testing/7814/7814gascomp.dat')
g = np.asarray(gas['zz'])

#***************************define total curve
#D=9.25 #disk M-L ratio provided in [1] 
#B=.5 #bulge M-L ratio provided in [1] 



# In[3]:


rcut=2.5 #default value
rho00=1e8 #default value

v_err1=v_err1
weighdata=1/v_err1
# LMFit
B=4.98331928
D=3.85244293
G=1.00239573
b=b*B
d=d*D
g=g*G
#r=np.linspace(minkpc,maxkpc,int(Max)) #kpc
