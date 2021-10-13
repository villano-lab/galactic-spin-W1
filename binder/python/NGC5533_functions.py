################################
########### Imports ############
################################
import sys
import traceback
import numpy as np
import scipy.special as ss
import scipy.optimize as so
import scipy.integrate as si
import scipy.interpolate as inter
try:
    import h5py as h5
    h5py = 1
except ModuleNotFoundError:
    h5py = 0
    print("Could not find h5py. Datasets will not be able to be saved or loaded using NGC5533_functions.")
#-----------For path detection-----------
import subprocess
#def getGitRoot():
#    return subprocess.Popen(['git', 'rev-parse', '--show-toplevel'], stdout=subprocess.PIPE).communicate()[0].rstrip().decode('utf-8')
# The above does not work on Windows machines and has therefore been removed.
defaultpath = '../'

################################
########## Constants ###########
################################

#---------Definitely Constant---------
G = 4.30091e-6    #gravitational constant (kpc/solar mass*(km/s)^2)
rhocrit = 9.3e-18 #critical density of the Universe (kg/km^3)

#---------Measured Directly-----------
L = 3.27e10                              #luminosity (Solar Luminosities)
absmag = -22.02                          #absolute magnitude
magsun = 4.42                            #absolute magnitude of the sun
L0 = np.power(10, (0.4*(magsun-absmag))) #Absolute Magnitude to luminosity

#---------Measured Indirectly---------
ups = 2.8                         #bulge mass-to-light ratio (Solar Mass/Solar Luminosity)???
q = 0.33                          #intrinsic axis ratio
e2 = 1-(q**2)                     #eccentricity
i = 52*(np.pi/180)                #inclination angle
h_rc = 1.4                        #core radius (kpc)
c = 1e-12                         #(what does this constant do?)
Mvir = 1e11*((c/(11.7))**(-40/3)) #virial mass (in solar mass) solved from eq(5)
Mbh_def = 2.7e9                   #Black Hole mass (in solar mass)

#---------Definitely Variable---------
n_c = 2.7                         #concentration parameter
h_c = 8.9                         #radial scale-length (kpc)
hrho00_c = 0.31e9                 #halo central surface density (solar mass/kpc^2)
drho00_c = 0.31e9                 #disk central surface density (solar mass/kpc^2)

#---------Uncategorized-------------------
re_c = 2.6                                               #effective radius (kpc)
epsdisk = 5.0                                            #from Noordermeer's paper
rs = (1/c)*(((3*Mvir)/((4*np.pi*100*rhocrit)))**(1/3))   #scale radius (kpc)
rho_s = (100/3)*((c**3)/(np.log(1+c)-(c/(1+c))))*rhocrit #characteristic density
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

################################
######### Black Hole ###########
################################

def bh_v(r,M=Mbh_def,save=False,load=True,comp='blackhole',**kwargs): #M in solar masses, r in kpc
    if isinstance(r,float) or isinstance(r,int):
        r = np.asarray([r])
    if isinstance(r,list):
        r = np.asarray(r)
    a = np.sqrt(G*M/r)
    if save:
        load = False
    if load:
        try: #Load existing prefactor if available
            #file = comp+'.hdf5'
            #print(file) #Troubleshooting
            y = loaddata(comp,'Mbh'+str(M),file=comp+'.hdf5',**kwargs)[1]
            x = loaddata(comp,'Mbh'+str(M),file=comp+'.hdf5',**kwargs)[0]
        except KeyError: #If unable to load, save
	   #If unable to load, load default instead and apply a prefactor retroactively
           # y = np.sqrt(M)*loaddata(comp,'Mbh1',file=comp+'.hdf5',**kwargs)[1]
           # x = loaddata(comp,'Mbh1',file=comp+'.hdf5',**kwargs)[0]
           # spline = inter.InterpolatedUnivariateSpline(x,y,k=3) #k is the order of the polynomial
           # return spline(r)
	        save = True
       # except: #Attempting to catch problem with spline having too few points
           # print('An error has occured. Switching to save function. Error information below:')
           # print(sys.exc_info()[0])
           # print(sys.exc_info()[1])
           # print()
           # print('#--------------------')
           # print()
           # print()
           # print(traceback.format_exc())
           # print()
           # print()
           # print('#--------------------')
           # print()
           # save = True #Calculate since there aren't enough points
    if save:
        savedata(r,a,comp,'Mbh'+str(M),file=comp+'.hdf5',**kwargs)
        return a
    else:
        return a
    
################################
########### Bulge ##############
################################

#I'm not sure how many of these we need to be defined -- I kept everything that was called outside of another function.
#We can condense the number of functions once we know for certain if there are things we don't need again.

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

def b_v(r,n=n_c,re=re_c,save=False,load=True,comp='bulge',**kwargs):
    if isinstance(r,float) or isinstance(r,int):
        r = np.asarray([r])
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
    a = b_vsquarev(r,n,re)**(1/2)
    a[np.isnan(a)] = 0
    if save:
        savedata(r,a,comp,'n'+str(n)+'re'+str(re),file=comp+'.hdf5',**kwargs)
        return a
    else:
        return a

################################
############ Halo ##############
################################

def h_rhat(r,z):               #r-hat from Casertano eq(9)
    return np.sqrt((r**2)+(z**2))

def h_rho(r,rho00=hrho00_c,rc=h_rc): #Isothermal Density Profile
    return rho00*((1+((r/rc)**2))**(-1))
    
def h_vcasertano(r,z,rc=h_rc,rho00=hrho00_c,gamma=h_gamma):                       #Velocity
    v0h = lambda r,rho00,rc,z: np.sqrt(h_rho(r,rho00,rc)*4*np.pi*G*(h_rhat(r,z)**2)) #eq 9 casertano
    return v0h(r,rho00,rc,z)*((r/rc)**gamma)                                     #eq 10 casertano

def h_vjimenez(r,rc=h_rc,rho00=hrho00_c):
    return np.sqrt(4*np.pi*G*rho00*(rc**2)*(1-((rc/r)*np.arctan(r/rc))))

def h_vNFW(r,save=True,comp='hNFW',**kwargs):
    rho = lambda r: rho_s/((r/rs)*((1+r/rs)**2))
    f = lambda R: 4*np.pi*rho(R)*(R**2)          #NFW Density Profile
    mdm = lambda r: si.quad(f, 0, r)[0]          #M(r)
    vdm2 = lambda r: (G*mdm(r))/r                #v^2: GM(r)/r
    vdm2v = np.vectorize(vdm2)
    a = np.sqrt(vdm2v(r))
    a[np.isnan(a)] = 0
    if save:
        savedata(r,a,comp,'n'+str('PLACEHOLDER'),file=comp+'.hdf5',**kwargs)
        return a
    elif load:
        return loaddata(comp,'n'+str('PLACEHOLDER'),file=comp+'.hdf5',**kwargs)
    else:
        return a(r)

def h_viso(r,rc=h_rc,rho00=hrho00_c,load=True,save=False,comp='halo',**kwargs):   #h_v iso
    if isinstance(r,float) or isinstance(r,int):
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

h_v = h_viso

        
################################
############ Disk ##############
################################

#----- To Fit For --------
#h, rho00

#----- Multiples of h ----
z0 = lambda h: 0.2*h                        #half-thickness (kpc)
R = lambda h: 4*h                           #cut-off radius (kpc)
d = lambda h: 0.2*h                         #cut-off length upper limits (kpc)

#----- Functions ---------

def d_px(r,u,xi):       #Initial Function
    #Matches Casertano
    x = lambda r,u,xi: (r**2+u**2+xi**2)/(2*r*u)
    try:
        return x(r,u,xi)-np.sqrt(x(r,u,xi)**2-1)
    except ZeroDivisionError: #If dividing by zero, return infinity instead of error. (Mostly at 0)
        return np.nan #This will allow nan handling later

def d_rho0(r, h=h_c, d_rho00=drho00_c): #density piecewise function
    #Matches Casertano
    conditions = [r <= R(h),
                  (r > R(h)) & (r <= R(h)+d(h)), 
                  r > R(h)+d(h)]
    functions = [lambda r,h,d_rho00: d_rho00*np.exp(-r/h),
                 lambda r,h,d_rho00: d_rho00*np.exp(-R(h)/h)*(1-(r-R(h))/d(h)),
                 lambda r,h,d_rho00: 0]
    return np.piecewise(r, conditions, functions, h, d_rho00)

def d_durho0(r, h=h_c, d_rho00=drho00_c): #partial derivative of rho(u,xi)
    #Doesn't seem to be written explicitly, but should match Casertano (derivative should be accurate at least where the problem is at the cusp (where it is 0))
    #Casertano takes the derivative of rho(u,z) with respect to u; here we take the derivative of rho(r,z) with repsect to r, which is equivalent.
    conditions = [r <= R(h),
                  (r > R(h)) & (r <= (R(h)+d(h))),
                  r > (R(h)+d(h))]
    functions = [lambda r,h,d_rho00: -(1/h)*d_rho00*np.exp(-r/h),
                 lambda r,h,d_rho00: -(1/d(h))*d_rho00*np.exp(-R(h)/h),
                 lambda r,h,d_rho00: 0]
    return np.piecewise(r, conditions, functions, h, d_rho00)

#Disk Density Distribution
def d_rho_rz(r,z,h=h_c,d_rho00=drho00_c):
    #Matches Casertano
    return d_rho0(r, h, d_rho00)*np.power(np.cosh(z/z0(h)), -2)
def d_drho_rz(r,z,h=h_c,d_rho00=drho00_c):
    #Doesn't seem to be written explicitly, but should match Casertano (derivative with respect to u, cosh term has no "u")
    try:
        return d_durho0(r, h, d_rho00)*np.power(np.cosh(z/z0(h)), -2)
    except ZeroDivisionError:
        return nan

def d_K(r,u,xi): #Complete Elliptic Integral
    #Matches Casertano
    return 2*(ss.ellipk(d_px(r,u,xi)) - ss.ellipe(d_px(r,u,xi)))/(np.pi*np.sqrt(r*u*d_px(r,u,xi)))

def d_innerfunc(z,r,u,h=h_c,d_rho00=drho00_c):  #Inner Function (3D)
    #Matches Casertano, where K = 2*(stuff)
    return d_drho_rz(u, z, h, d_rho00)*d_K(r,u,z)

def d_innerintegral(u,r,h=h_c,d_rho00=drho00_c): #Integrate Function
    #Matches Casertano, with the exception of the 125 limit due to computing limitations
    #u was passed from outerintegral to prevent errors.
    return u*si.quad(d_innerfunc, 0, R(h)+d(h), args=(r,u,h,d_rho00))[0]
#Args passed into quad need to be numbers, not arrays. (?)

def d_outerintegral(r,h=h_c,d_rho00=drho00_c): #Integrate Outer Function
    #Matches casertano. u was passed at the end of innerintegral to prevent errors.
    return si.quad(d_innerintegral, 0, np.inf, args=(r,h,d_rho00))[0]

def d_Mdblintrho(h=h_c,d_rho00=drho00_c):    #M double-integral rho
    rho_rz_r = lambda z,r,h,d_rho00: d_rho_rz(r,z,h,d_rho00)*r
    Mintrho = lambda r,h,d_rho00: si.quad(rho_rz_r, -(R(h)+d(h)), R(h)+d(h), args=(r,h,d_rho00))[0]
    return si.quad(Mintrho,0,np.inf,args=(h,d_rho00))[0]

pref_def = 2.352579926191481 #epsdisk*(L0/d_Mdblintrho(defaults))

def d_F(r,h=h_c,d_rho00=drho00_c,pref=1): #multiplying by upsylon. Generally we should either use pref or h and d
    #Matches Casertano
    if pref == False:
        pref = epsdisk*(L0/d_Mdblintrho(h,d_rho00))
    val = 4*np.pi*G*d_outerintegral(r,h,d_rho00)*pref
    if np.isnan(val):
        return 0
    else:
        return val
d_Fv = np.vectorize(d_F)

#NOTE: The order of variables for d_v is different than above, and is different from how they are listed in the load/save.
#This should not adversely affect anything, just be careful that you have inputs in the correct order for the given function.
#This was done for conveneince; pref is now first and easy to call as the only non-default.
def d_v(r,pref=0.5,h=h_c,d_rho00=drho00_c,save=False,load=True,comp='disk',**kwargs): #velocity
    #Matches Casertano
    if isinstance(r,float) or isinstance(r,int):
        r = np.asarray([r])
    if isinstance(r,list):
        r = np.asarray(r)
    if load:
        try: #Load existing prefactor if available
            y = loaddata(comp,'h'+str(h)+'d_rho00'+str(d_rho00)+'pref'+str(pref),file=comp+'.hdf5',**kwargs)[1]
            x = loaddata(comp,'h'+str(h)+'d_rho00'+str(d_rho00)+'pref'+str(pref),file=comp+'.hdf5',**kwargs)[0]
            b = inter.InterpolatedUnivariateSpline(x,y,k=3) #k is the order of the polynomial
            return b(r)
        except KeyError: #If unable to load, load 1 instead and apply a prefactor retroactively
            try:
                y = pref*loaddata(comp,'h'+str(h)+'d_rho00'+str(d_rho00)+'pref1',file=comp+'.hdf5',**kwargs)[1]
                x = loaddata(comp,'h'+str(h)+'d_rho00'+str(d_rho00)+'pref1',file=comp+'.hdf5',**kwargs)[0]
                b = inter.InterpolatedUnivariateSpline(x,y,k=3) #k is the order of the polynomial
                return b(r)
            except KeyError: #And if still unable to load, calculate and save.
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
        r = np.asarray(r)
        a = np.sqrt(-r*d_Fv(r,h,d_rho00,pref))
        a[np.isnan(a)] = 0
        savedata(r,a,comp,'h'+str(h)+'d_rho00'+str(d_rho00)+'pref'+str(pref),file=comp+'.hdf5',**kwargs)
        return a
    else:
        a = np.sqrt(-r*d_Fv(r,h,d_rho00,pref))
        a[np.isnan(a)] = 0
        return a

def d_thief(r,pref=1):
    sys.path.append('./')
    import dataPython as dp
    data = dp.getXYdata('../data/final/nord-120kpc-disk.txt')
    x = np.asarray(data['xx'])
    y = np.asarray(data['yy'])
    b = inter.InterpolatedUnivariateSpline(x,y,k=1) #k is the order of the polynomial    
    return b(r)

def g_thief(r,pref=1):
    sys.path.append('./')
    import dataPython as dp
    data = dp.getXYdata('../data/final/nord-120kpc-gas.txt')
    x = np.asarray(data['xx'])
    y = np.asarray(data['yy'])
    b = inter.InterpolatedUnivariateSpline(x,y,k=1) #k is the order of the polynomial    
    return b(r)
    
   
    
    
################################
############ Total #############
################################
def v(r,M=Mbh_def,re=re_c,h=h_c,d_rho00=drho00_c,pref=1,rc=h_rc,h_rho00=hrho00_c,save=False,load=True,comp='total',**kwargs): 
    if isinstance(r,float) or isinstance(r,int):
        r = np.asarray([r])
    if load:
        try: #Load if exists
            y = loaddata(comp,'Mbh'+str(M)+'re'+str(re)+'h'+str(h)+'d_rho00'+str(d_rho00)+'pref'+str(pref) +'rc'+str(rc)+'h_rho00'+str(h_rho00), file=comp+'.hdf5',**kwargs)[1]
            x = loaddata(comp,'Mbh'+str(M)+'re'+str(re)+'h'+str(h)+'d_rho00'+str(d_rho00)+'pref'+str(pref) +'rc'+str(rc)+'h_rho00'+str(h_rho00), file=comp+'.hdf5', **kwargs)[0]
            b = inter.InterpolatedUnivariateSpline(x,y,k=3)
            return b(r)
        except KeyError: #If does not exist,
            save = True  #Save instead
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
    a = np.sqrt(np.sqrt(h_v(r,rc,h_rho00,load,save)**2+d_v(r,h,d_rho00,pref,load,save)**2+bh_v(r,M,load,save)**2+b_v(r,re,load,save)**2))
    a[np.isnan(a)] = 0
    if save: #not elif since that would mean don't check if load was true, which I don't want in this case
        savedata(r,a,comp,'Mbh'+str(M)+'re'+str(re)+'h'+str(h)+'d_rho00'+str(d_rho00)+'pref'+str(pref) +'rc'+str(rc)+'h_rho00'+str(h_rho00),file=comp+'.hdf5',**kwargs)
        return a
    elif not load: #If load was false,
        return a
