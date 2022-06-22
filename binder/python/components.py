###############
### Imports ###
###############

import numpy as np
import dataPython as dp
import inspect
import types
import sys
import scipy.integrate as si
import scipy.optimize as so
import scipy.special as ss
import scipy.interpolate as inter
from scipy.interpolate import InterpolatedUnivariateSpline      # Spline function
import lmfit as lm                                              # Fitting
#custom libraries
import load_galaxies as galdata

try:
    import h5py as h5
    h5py = 1
except ModuleNotFoundError:
    h5py = 0
    print("Could not find h5py. Datasets will not be able to be saved or loaded using components.py.")

defaultpath = '../'
ngc7814list = ['ngc7814', 'ngc 7814', '7814']
ngc5533list = ['ngc5533', 'ngc 5533', '5533']
ngc0891list = ['ngc0891', 'ngc 0891', '0891', 'ngc891', 'ngc 891', '891']

################################
########## Constants ###########
################################
#Pull from load_galaxies.py dict
def dict(galaxy):
    if galaxy.lower() in ngc7814list:
        return galdata.NGC7814
    if galaxy.lower() in ngc5533list:
        return galdata.NGC5533
    if galaxy.lower() in ngc0891list:
        return galdata.NGC0891

#Defaults based on NGC5533

#---------Definitely Constant---------
G = 4.30091e-6                           # Gravitational constant (kpc/solar mass*(km/s)^2)
#rhocrit = 9.3e-18                        # Critical density of the Universe (kg/km^3). 

#---------Measured Directly-----------
L = 3.27e10                              # Luminosity (Solar Luminosities)
absmag = -22.02                          # Absolute magnitude
magsun = 4.42                            # Absolute magnitude of the sun
L0 = np.power(10, (0.4*(magsun-absmag))) # Absolute magnitude to luminosity

#---------Measured Indirectly---------
ups = 2.8                         # Bulge mass-to-light ratio (Solar Mass/Solar Luminosity). Source: Noordermeer, 2008
q = 0.33                          # Intrinsic axis ratio. Source: Noordermeer, 2008
e2 = 1-(q**2)                     # Eccentricity. Source: Noordermeer, 2008
i = 52*(np.pi/180)                # Inclination angle. Source: Noordermeer & Van Der Hulst, 2007
h_rc = 1.4                        # Core radius (kpc). Source: Noordermeer, 2008
#c = 1e-12                         #(what does this constant do?)
#Mvir = 1e11*((c/(11.7))**(-40/3)) # Virial mass (in solar mass) solved from eq(5)
Mbh_def = 2.7e9                   # Black Hole mass (in solar mass). Source: Noordermeer, 2008

#---------Definitely Variable---------
n_c = 2.7                         # Concentration parameter. Source: Noordermeer & Van Der Hulst, 2007
h_c = 8.9                         # Radial scale-length (kpc). Source: Noordermeer & Van Der Hulst, 2007
hrho00_c = 0.31e9                 # Halo central surface density (solar mass/kpc^2)
drho00_c = 0.31e9                 # Disk central surface density (solar mass/kpc^2)

#---------Uncategorized-------------------
re_c = 2.6                        # Effective radius (kpc). Source: Noordermeer & Van Der Hulst, 2007
upsdisk = 5.0                     # Disk mass-to-light ratio. Source: Noordermeer, 2008
#rs = (1/c)*(((3*Mvir)/((4*np.pi*100*rhocrit)))**(1/3))   # Scale radius (kpc)
#rho_s = (100/3)*((c**3)/(np.log(1+c)-(c/(1+c))))*rhocrit # Characteristic density
h_gamma = 0

################################
########### Saving #############
################################

def savedata(xvalues,yvalues,group,dataset,path=defaultpath,file='Inputs.hdf5'): 
#this is a dummy filename to enforce ordering; try not to save here except for testing!
    if h5py == 1:
        saved = h5.File(path+'/'+file,'a')
        if group in ['Disk', 'disc', 'Disc',  'd', 'D']:
            group = 'disk'
            print("Group name set to 'disk'.")
        if group in ['bh','Bh','BH','Black Hole','BlackHole','Blackhole,','Black hole','black hole','Black Hole']:
            group = 'blackhole'
            print("Group name set to 'blackhole'.")
        if group in ['dm','DM','Dm','Dark Matter','Dark matter','dark matter','h','H','Halo','darkmatter','Darkmatter','DarkMatter']:
            group = 'halo'
            print("Group name set to 'halo'.")
        if group in ['b','B','Bulge']:
            group = 'bulge'
            print("Group name set to 'bulge'.")
        if group in ['t','T','Total']:
            group = 'total'
            print("Group name set to 'total'.")
        try:
            grp = saved.create_group(group)
            grp.create_dataset(dataset,data=[xvalues,yvalues])
        except ValueError:
            try:
                grp = saved[group]
                grp.create_dataset(dataset,data=[xvalues,yvalues])
            except RuntimeError:
                x = loaddata(group,dataset,path,file)[0]
                x = np.append(x,xvalues)
                y = loaddata(group,dataset,path,file)[1]
                y = np.append(y,yvalues)
                x, y = (list(a) for a in zip(*sorted(zip(x, y))))
                i = 0
                while i < len(x)-1:
                    if x[i+1] == x[i]:
                        x = np.delete(x,i+1)
                        y = np.delete(y,i+1)
                    else:
                        i += 1
                del grp[dataset]
                savedata(x,y,group,dataset,path,file)
                return y
        finally: #No matter what,
            saved.close()
        #print("Saved.") #Convenient for debugging but annoying for fitting.
    elif h5py == 0:
        print("ERROR: h5py was not loaded.")
        return 1
    
def loaddata(group,dataset,path=defaultpath,file='Inputs.hdf5'):
    if h5py == 1:
        saved = h5.File(path+'/'+file,'r')
        if group in ['Disk', 'disc', 'Disc', 'd', 'D']:
            group = 'disk'
            print("Group name set to 'disk'.")
        if group in ['bh','Bh','BH','Black Hole','BlackHole','Blackhole,','Black hole','black hole','Black Hole']:
            group = 'blackhole'
            print("Group name set to 'blackhole'.")
        if group in ['dm','DM','Dm','Dark Matter','Dark matter','dark matter','h','H','Halo','darkmatter','Darkmatter','Dark Matter']:
            group = 'halo'
            print("Group name set to 'halo'.")
        if group in ['b','B','Bulge']:
            group = 'bulge'
            print("Group name set to 'bulge'.")
        if group in ['t','T','Total']:
            group = 'total'
            print("Group name set to 'total'.")
        grp = saved[group]
        dset = grp[dataset]
        a = dset[:]
        return a
    #Placeholder; I will design this to store information at a later date.
    elif h5py == 0:
        print("ERROR: h5py was not loaded.")
        return 1
    saved.close() #no matter what, close the file when you're done
    
def checkfile(group='all',path=defaultpath,file='Inputs.hdf5'):
    if h5py ==1:
        saved = h5.File(path+'/'+file,'r')
        if group == 'all':
            print('Groups:')
            for n in saved:
                print(saved[n])
            print('')
            print(' ---------------- ')
            print('')
            print('More information:')
            for n in saved:
                grp = saved[n]
                print(grp)
                for m in grp:
                    print('        '+str(grp[m]))
        else:
            print(group+':')
            grp = saved[group]
            for n in grp:
                print(grp[n])
        saved.close()
    elif h5py ==0:
        print("ERROR: h5py was not loaded.")
        return 1

    
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
        return galdata.data
def data_total(galaxy):
    if galaxy.lower() in ngc5533list:
        return galdata.data_total
    else:
        return None
def r_dat(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(data(galaxy)['xx'])
    elif galaxy.lower() in ngc5533list:
        return galdata.NGC5533['m_radii']
def v_dat(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(data(galaxy)['yy'])
    elif galaxy.lower() in ngc5533list:
        return galdata.NGC5533['m_velocities']
def v_err0(galaxy):
    if galaxy.lower() in ngc5533list:
        return galdata.v_err0
    else:
        return None
def v_err1(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(data(galaxy)['ey'])
    elif galaxy.lower() in ngc5533list:
        return galdata.NGC5533['m_v_errors']
    else:
        raise ValueError('Galaxy string', galaxy, 'not recognized.')
    
def disk_raw(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(dp.getXYdata('data/NGC7814/7814_gDisk.dat')['yy'])
    else:
        return None

def bulge_raw(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(dp.getXYdata('data/NGC7814/7814_gBulge.dat')['yy'])
    else:
        return None
    
def gas_raw(galaxy):
    if galaxy.lower() in ngc7814list:
        return np.asarray(dp.getXYdata('data/NGC7814/7814_gGas.dat')['yy'])
    else:
        return None
    
#####################
### Interpolation ###
#####################

def interpd(x,y):
    return InterpolatedUnivariateSpline(x,y,k=5)

################################
######### Components ###########
################################

def blackhole(r,M,galaxy,load=False,save=False):
    if galaxy.lower() in ngc5533list:
        if isinstance(r,float) or isinstance(r,int):
            r = np.asarray([r])
        if isinstance(r,list):
            r = np.asarray(r)
        a = np.sqrt(G*M/r)
        if save:
            load = False
        if load:
            try: #Load existing prefactor if available
                y = loaddata(comp,'Mbh'+str(M),file=comp+'.hdf5',**kwargs)[1]
                x = loaddata(comp,'Mbh'+str(M),file=comp+'.hdf5',**kwargs)[0]
            except KeyError: #If unable to load, save
                save = True
        if save:
            savedata(r,a,comp,'Mbh'+str(M),file=comp+'.hdf5',**kwargs)
            return a
        else:
            return a
    else: #galaxy has no significant bh component
        return 0
    
def b_gammafunc(x,n=n_c):
    return ss.gammainc(2*n,x)*ss.gamma(2*n)-0.5*ss.gamma(2*n)
b_root = so.brentq(b_gammafunc,0,500000,rtol=0.000001,maxiter=100) #come within 1% of exact root within 100 iterations
def b_I0(n=n_c,re=re_c):
    return L*(b_root**(2*n))/(re**2*2*np.pi*n*ss.gamma(2*n))
def b_r0(n=n_c,re=re_c):
    return re/np.power(b_root,n)
def b_innerintegral(m,n=n_c,re=re_c):
    f = lambda x,m,n,re: np.exp(-np.power(x/b_r0(n,re), (1/n)))*np.power(x/b_r0(n,re), 1/n-1)/(np.sqrt(x**2-m**2)) #Inner function
    return si.quad(f, m, np.inf,args=(m,n,re))[0]
b_innerintegralv = np.vectorize(b_innerintegral)
def b_vsquare(r,n=n_c,re=re_c):
    C = lambda n,re: (4*G*q*ups*b_I0(n,re))/(b_r0(n,re)*np.float(n))*(np.sqrt((np.sin(i)**2)+(1/(q**2))*(np.cos(i)**2)))
    h = lambda m,r,n,re: C(n,re)*b_innerintegral(m,n,re)*(m**2)/(np.sqrt((r**2)-((m**2)*(e2)))) #integrate outer function
    return si.quad(h, 0, r, args=(r,n,re))[0]
def b_vsquarev(r,n=n_c,re=re_c):
    a = np.vectorize(b_vsquare)
    return a(r,n,re)
def bulge(r,bpref,galaxy,n=n_c,re=re_c,load=True,save=False,**kwargs):
    if galaxy.lower() in ngc7814list:
        polynomial = interpd(r_dat(galaxy),bpref*bulge_raw(galaxy))   
        return polynomial(r)
    elif galaxy.lower() in ngc5533list:
        x = np.sort(r_dat(galaxy))
        if isinstance(r,float) or isinstance(r,int):
            r = np.asarray([r])
        comp = 'bulge'
        if load:
            try: #load if exists
                y = loaddata(comp,'n'+str(n)+'re'+str(re),file=comp+'.hdf5',**kwargs)[1]
                x = loaddata(comp,'n'+str(n)+'re'+str(re),file=comp+'.hdf5',**kwargs)[0]
                b = inter.InterpolatedUnivariateSpline(x,y,k=3) #k is the order of the polynomial
                return b(r)
            except KeyError: #if does not exist,
                save = True  #go to save function instead
            except: #Attempting to catch problem with spline having too few points
                print('An error has occured. Switching to save function. Error information below:')
                print(sys.exc_info()[0])
                print(sys.exc_info()[1])
                print()
                print('#--------------------')
                print()
                print()
                print(traceback.format_exc())
                print()
                print()
                print('#--------------------')
                print()
                save = True #Calculate since there aren't enough points
        a = b_vsquarev(r_dat(galaxy),n,re)**(1/2)
        a[np.isnan(a)] = 0
        if save:
            savedata(r,a,comp,'n'+str(n)+'re'+str(re),file=comp+'.hdf5',**kwargs)
        y = bpref*a
        polynomial = interpd(x,y)
        return polynomial(r)

def disk(r,dpref,galaxy):
    if galaxy.lower() in ngc7814list:
        polynomial = interpd(r_dat(galaxy),dpref*disk_raw(galaxy))   
        return polynomial(r)
    elif galaxy.lower() in ngc5533list:
        data = dp.getXYdata('data/NGC5533/noord-120kpc-disk.txt')
        x = np.asarray(data['xx'])
        y = np.asarray(data['yy'])
        d = inter.InterpolatedUnivariateSpline(x,y,k=3) #k is the order of the polynomial
        return dpref*d(r)

def gas(r,gpref,galaxy):
    if galaxy.lower() in ngc7814list:
        polynomial = interpd(r_dat(galaxy),gpref*gas_raw(galaxy))   
        return polynomial(r)
    elif galaxy.lower() in ngc5533list:
        data = dp.getXYdata('data/NGC5533/noord-120kpc-gas.txt')
        x = np.asarray(data['xx'])
        y = np.asarray(data['yy'])
        g = inter.InterpolatedUnivariateSpline(x,y,k=3) #k is the order of the polynomial    
        return gpref*g(r)
    
#########################
### Galaxy parameters ###
#########################


def rcut(galaxy):
    if galaxy.lower() in ngc7814list: # NGC 7814 (Source: Fraternali, Sancisi, and Kamphuis, 2011)
        return 2.1            # Cutoff (core) radius (in kpc) from Table 5.
    elif galaxy.lower() in ngc5533list:
        return 1.4 # Core radius (kpc). Source: Noordermeer, 2008
    
def rho0(galaxy):
    if galaxy.lower() in ngc7814list:
        return 1.524e8            # Central density (in solar mass/kpc^3) from Table 5. (converted from pc to kpc)
    elif galaxy.lower() in ngc5533list:
        return 0.31e9                 # Halo central surface density (solar mass/kpc^2)

# Constants
G = 4.30091e-6            # Gravitational constant (kpc/solar mass*(km/s)^2)

#################################
### Calculating enclosed mass ### 
#################################

# Mass as a function of radius, calculated from the Isothermal density profile:
def mass_r(r,rcut):
    return 4 * np.pi * rcut**2  * (r - rcut * (np.arctan(r/rcut))) 
# rho0 represents the number of tiny black holes at the center of the galaxy

######################################################
### Calculating halo velocity with various methods ###
######################################################

# Calculating the velocity for each black hole as a point mass
def halo_BH(r,scale,arraysize,massMiniBH,rcut):
    x = np.sort(r)
    y = np.sqrt((G * (scale * arraysize * massMiniBH) * mass_r(r,rcut)) / r)   
                    # scale is needed to be separate and constant because the widget would freeze the computer otherwise
                    # arraysize is the number of black holes for slider
                    # mBH is the mass of black holes for slider
    halo = interpd(x,y)
    return halo(r)

def h_viso(r,rc=rcut('5533'),rho00=rho0('5533'),load=True,save=False,comp='halo',**kwargs):   #h_v iso
    if isinstance(r,float) or isinstance(r,int): #if r isn't array-like, make it array-like.
        r = np.asarray([r])
    a = np.zeros(len(r))
    i = 1
    while i < len(r):
        a[i] = np.sqrt(
            4*np.pi*G*rho00*(rc**2)*(1-(
                (rc/r[i])*np.arctan(r[i]/rc))
                                    )
        )
        i += 1
    a[np.isnan(a)] = 0
    if load:
        try: #Load if exists
            y = loaddata(comp,'rc'+str(rc)+'rho00'+str(rho00),file=comp+'.hdf5',**kwargs)[1]
            x = loaddata(comp,'rc'+str(rc)+'rho00'+str(rho00),file=comp+'.hdf5',**kwargs)[0]
            b = inter.InterpolatedUnivariateSpline(x,y,k=3) #k is the order of the polynomial
            return b(r)
        except KeyError: #If does not exist,
            save = True #Calculate and save
        except FileNotFoundError:
            save = True
        except: #Attempting to catch problem with spline having too few points
            print('An error has occured. Switching to save function. Error information below:')
            print(sys.exc_info()[0])
            print(sys.exc_info()[1])
            print()
            print('#--------------------')
            print()
            print()
            print(traceback.format_exc())
            print()
            print()
            print('#--------------------')
            print()
            save = True #Calculate since there aren't enough points
    if save:
        savedata(r,a,comp,'rc'+str(rc)+'rho00'+str(rho00),file=comp+'.hdf5',**kwargs)
        return a
    else:
        return a

def halo(r,rc,rho00): #A 'default' version
        return h_viso(r,rc,rho00,load=False)

##################################
### Calculating total velocity ###
##################################
def totalvelocity_miniBH(r,scale,arraysize,massMiniBH,rcut,bpref,dpref,gpref,Mbh,galaxy):
    return np.sqrt(blackhole(r,Mbh,galaxy)**2 
                        + bulge(r,bpref,galaxy)**2 
                        + disk(r,dpref,galaxy)**2
                        + halo_BH(r,scale,arraysize,massMiniBH,rcut)**2
                        + gas(r,gpref,galaxy)**2)
    
def totalvelocity_halo(r,scale,arraysize,rho00,rcut,bpref,dpref,gpref,Mbh,galaxy):
    return np.sqrt(blackhole(r,Mbh,galaxy,load=False)**2
                   + bulge(r,bpref,galaxy)**2
                   + disk(r,dpref,galaxy)**2
                   + halo(r,rcut,rho00)**2
                   + gas(r,gpref,galaxy)**2)
    
#################################
#### Find Fitting Parameters ####
#################################
# For tiny black hole widget

# Setup
def weighdata(galaxy):
    return 1/v_err1(galaxy)

# Tiny Black Holes 
def set_params(model,galaxy):
    if type(model) == types.FunctionType:
        model = lm.Model(model)
        fit_pars = model.make_params()
    elif type(model) == lm.Model:
        fit_pars = model.make_params()
    else:
        raise ValueError("Invalid type for variable `model`. (",model,", type:",type(model),".)")
    #miniBH halo
    #fit_pars.add('rc',    value=rcut(galaxy), min=0.1)         # Core Radius (kpc)
    fit_pars.add('scale',      value=rho0(galaxy),   vary=False)        # Scale
    if (model == 'bh') or (model==lm.Model(totalvelocity_miniBH)) or (model==totalvelocity_miniBH):
        fit_pars.add('arraysize',  value=50,     min=1, max=100)    # Number of black holes
        fit_pars.add('rho0', value=1.5, min=0)       # Halo Central Density 
    elif (model == 'wimp') or (model==lm.Model(totalvelocity_halo)) or (model==totalvelocity_halo):
        fit_pars.add('arraysize', value=0, vary=False)
        fit_pars.add('rho0', value=rho0(galaxy), min=0)
    fit_pars.add('rcut',       value=rcut(galaxy),   min=0.1)           # Core Radius (kpc)
    # Bulge
    fit_pars.add('bpref', value=1, min=0, max=100)  # Bulge Prefactor
    # Disk
    fit_pars.add('dpref', value=1, min=0, max=100)  # Disk Prefactor
    # Disk
    fit_pars.add('gpref', value=1, vary=False)        # Gas Prefactor
    # BH
    if galaxy in ngc5533list:
        fit_pars.add('Mbh', value=2.7e9, min=1e8) # Black Hole mass (in solar mass). Source: Noordermeer, 2008
    else:
        fit_pars.add('Mbh', value=0, vary=False)
    return fit_pars

# Do fit
def bestfit(model,galaxy):
    newmodel = lambda r,scale,arraysize,rho0,rcut,bpref,dpref,gpref,Mbh: model(r,scale,arraysize,rho0,rcut,bpref,dpref,gpref,Mbh,galaxy)
    fit_mod = lm.Model(newmodel)
    fit_pars = set_params(newmodel,galaxy)
    if (model == 'bh') or (model==lm.Model(totalvelocity_miniBH)) or (model==totalvelocity_miniBH):
        fit_pars.add('arraysize',  value=50,     min=1, max=100)    # Number of black holes
        fit_pars.add('rho0', value=1.5, min=0)       # Halo Central Density 
    elif (model == 'wimp') or (model==lm.Model(totalvelocity_halo)) or (model==totalvelocity_halo):
        fit_pars.add('arraysize', value=0, vary=False)
        fit_pars.add('rho0', value=rho0(galaxy), min=0)
    if galaxy.lower() in ngc5533list and (model == lm.Model(totalvelocity_halo) or model == totalvelocity_halo):
        weights = 1/np.sqrt(v_err1(galaxy)**2+galdata.NGC5533['n_v_bandwidth']**2)
    else:
        weights = weighdata(galaxy)
    fit = fit_mod.fit(v_dat(galaxy),fit_pars,r=r_dat(galaxy),weights=weights)
    bestfit = fit.best_fit
    fit_dict = fit.best_values
    return bestfit, fit_dict
    #Examples on how to use this output:
    #best_rc = fit_dict['rc']
    #best_rho00 = fit_dict['rho00']
    #best_bpref = fit_dict['bpref']
    #best_dpref = fit_dict['dpref']
    #best_gpref = fit_dict['gpref']