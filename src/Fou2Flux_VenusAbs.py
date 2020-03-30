#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:48:57 2018

Script to calculate F and P from a set of fourier files 
and a given set of geometry. 

this script also has plotting functions for variation of band depth with wavelength. 

@author: gouravmahapatr, TU Delft
"""
import os
import pymiedap.pymiedap as pmd
import pymiedap.data_class as dt
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import matplotlib.style
import matplotlib as mpl
import linestyles as ls
mpl.style.use('classic') # set the default plotting style to 'classic'

# first load all the geometry information related to a particular orbit from the VEX data
path = '/Users/gouravmahapatr/Dropbox/PhD/spicav_data'
os.chdir(path)
import read_spicav_ir as rd

"""
load a specific orbit on the basis of its orbit number 
"""
data = rd.processing_orbits(1478,8)

def geom(index=-1):
    '''
    read the orbit geometries for a pass.
    Format: [phase,emission,sza,lat,lon]
    '''
    phase0 = data.geo.phase
    emission0 = data.geo.emission
    sza0 = data.geo.sza
    lat0 = data.geo.Lat
    lon0 = data.geo.Lon

    #phase = phase0[-2:-1]
    #emission = emission0[-2:-1]
    #sza = sza0[-2:-1]
    phase = phase0[np.array(index)]
    emission = emission0[np.array(index)]
    sza = sza0[np.array(index)]
    lat = lat0[np.array(index)]
    lon = lon0[np.array(index)]
    
    geo = [np.array([phase]),np.array([emission]),np.array([sza]),np.array([lat]),np.array([lon])]
    
    return np.array(geo)


# load the bmabs data used for generating the fourier files
#path = '/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/DAPlbl1/'
#path = '/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/DAPvenusCombLay'
#os.chdir(path)

#kdgw = np.loadtxt('kdgauss.out')

#xpDeep = np.linspace(-3,6,100)
xpDeep = np.linspace(-2,6,400)
absTDeep = 10**xpDeep  # this is the total bmabs array

def absDeep(lamb):
    bmabs = np.loadtxt('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/DAPabs1/bmabsDeep.dat',skiprows=4)
    #bmabs = np.loadtxt('bmabs.dat',skiprows=4)
    #
    #lambda1 = np.loadtxt('wav.dat')
    #xp1 = np.linspace(-3,4,200)
    #xp2 = np.linspace(3,6,200)
    #xpCont = np.zeros(371)
    #xpCont[0:171] = xp1[0:171]
    #xpCont[171:371] = xp2
    #absTCont = 10**xpCont
      
    # get the geometry for which the Stokes will be generated
    geo = geom()
    geo[0,0] = 30
    geo[2,0] = 30
    geo[1,0] = 0
    
    ## change the path to the location of fourier files
    path = '/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/DAPabs1/FourierFilesCloud+Haze0.2@{:04d}nm'.format(lamb)
    ##path = '/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/DAPlbl1/FourierFiles_lbl_1.44t1.45'
    os.chdir(path)
    #
    ## initialize the stokes arrays for multiple absorption values
    Iar = np.zeros((bmabs.shape[0],geo.shape[1]))
    Qar = np.zeros((bmabs.shape[0],geo.shape[1]))
    Uar = np.zeros((bmabs.shape[0],geo.shape[1]))
    Var = np.zeros((bmabs.shape[0],geo.shape[1]))
    Par = np.zeros((bmabs.shape[0],geo.shape[1]))
    
    # call the dap function from the pymiedap code and calculate the stokes vectors
    num = 1
    for i in range(bmabs.shape[0]-1):
        fname = 'fou'+str('%04d'%num)+'.dat'
        I,Q,U,V = pmd.read_dap_output(geo[0],geo[2],geo[1],fname,beta=[0],phi=[0])
        num+=1
        Iar[i,:] = np.double(I)
        Qar[i,:] = np.double(Q)
        Uar[i,:] = np.double(U)
        Var[i,:] = np.double(V)
        print(num)
    Iar_sm = savgol_filter(Iar[:,0],31,3)
    Qar_sm = savgol_filter(Qar[:,0],31,3)
    Uar_sm = savgol_filter(Uar[:,0],31,3)
    #Par = np.sqrt(Qar_sm**2 + Uar_sm**2)/Iar_sm
    Par = -Qar_sm/Iar_sm
    
    return Iar_sm,Qar_sm,Uar_sm,Par

# %% summon the function here for different wavelengths

I = np.zeros((xpDeep.shape[0],6))
Q = np.zeros((xpDeep.shape[0],6))
U = np.zeros((xpDeep.shape[0],6))
P = np.zeros((xpDeep.shape[0],6))

I[:,0],Q[:,0],U[:,0],P[:,0] = absDeep(870)
I[:,1],Q[:,1],U[:,1],P[:,1] = absDeep(950)
I[:,2],Q[:,2],U[:,2],P[:,2] = absDeep(1050)
I[:,3],Q[:,3],U[:,3],P[:,3] = absDeep(1200)
I[:,4],Q[:,4],U[:,4],P[:,4] = absDeep(1450)
I[:,5],Q[:,5],U[:,5],P[:,5] = absDeep(1750)

# %% set all the values of I smaller than eps to NaNs
eps = 1.e-5
#I[np.where(I<eps)] = np.nan
#P[np.where(I<eps)] = np.nan


# %% only for plotting

lamblist = ['870 nm','950 nm','1050 nm','1200 nm','1450 nm','1750 nm']
colors = ['purple','blue','green','orange','red','black']
plt.figure(figsize=[18,6])
lines = [ls.linestyles['solid'],ls.linestyles['dashed'],ls.linestyles['densely dashed'],ls.linestyles['dotted'],ls.linestyles['densely dotted'],ls.linestyles['dashdotted']]
findx = 300
for i in range(len(lamblist)):   
    plt.subplot(121)
    plt.semilogx(absTDeep[:findx],I[:findx,i],ls=lines[i],lw=1.5,c='k',label=lamblist[i])
    #if i == 5: plt.semilogx(absTDeep[:findx],I[:findx,i],ls=ls,lw=1.5,c='k')
    plt.subplot(122)
    plt.semilogx(absTDeep[:findx],P[:findx,i],ls=lines[i],lw=1.5,c='k')

plt.subplot(122);plt.ylim([-0.04,0.04]);plt.xlim([10**-2,10**3.2])
plt.subplot(121);plt.xlim([10**-2,10**3.2])
plt.subplot(121);plt.legend()
#plt.subplot(122);plt.legend(lamblist)
#plt.subplot(122);plt.legend(lamblist,loc=3)
plt.subplot(121)
plt.grid()
plt.ylabel('I',fontsize='large')
plt.xlabel('$\sum b^m_{abs}$',fontsize='large')
plt.subplot(122)
plt.grid()
plt.ylabel('P',fontsize='large')
plt.xlabel('$\sum b^m_{abs}$',fontsize='large')
plt.tight_layout()
    
    
# %% 
    
#plt.subplot(121)
#plt.xlim([10**-2,10**3])
#plt.legend(lamblist)
#plt.ylabel('I',fontsize='large')
#plt.grid()
#plt.subplot(122)
#plt.xlim([10**-2,10**3])
##plt.ylim([-0.05,0.05])
#plt.ylabel('P',fontsize='large')
#plt.grid()
#plt.subplot(223)
##plt.xlim([10**-2,10**3])
#plt.ylabel('Q',fontsize='large')
#plt.xlabel('$\sum b^m_{abs}$',fontsize='large')
#plt.grid()
#plt.subplot(224)
##plt.xlim([10**-2,10**3])
#plt.ylabel('U',fontsize='large')
#plt.xlabel('$\sum b^m_{abs}$',fontsize='large')
#plt.grid()
#plt.tight_layout()
#
##y = smooth(Qar[:,0],window_len=21)
##absT_sm = absTDeep
##Qar_sm = y
##Iar_sm = smooth(Iar[:,0],window_len=21)
#
#Par = -Qar[:,0]/Iar[:,0]
##Par = -Qar_sm/Iar[:,0]

def flarr():
    # enter the starting lambda, d(lambda) and number of lambda
    lamb0 = 1.40050
    dlamb = 0.0005
    nlamb = 200
    # make sure the number of decimal places for file reading is correct 
    larr_str = [str('%1.5f'%(lamb0+dlamb*i)) for i in range(nlamb)]
    larr = np.array(larr_str,dtype=float)
    
    return larr_str,larr
    
    
def foukd2flux():    
    # get the geometry of the observation, 
    # default geometry is closest to equator
    phase,emi,sza,lat,lon = geom(index=0)
    
    # point the function to the fourier folder location
    set_path0 = '/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/DAPvenusCombLay/FourierFiles'
    os.chdir(set_path0)
    
    larr,nlarr = flarr()
    
    I,Q,U,V = [],[],[],[]
    for j in range(len(larr)):
        print('Wavelength: ',larr[j])
        # enter the folder with the wavelength's name
        set_path1 = set_path0+'/'+larr[j]
        os.chdir(set_path1)
        
        Il,Ql,Ul,Vl = [],[],[],[]
        
        for i in range(10):
            fname = 'fou'+str(i+1)+'.out'
            #print('Now reading... ',fname)
            
            # use pymiedap's read_dap_output to get the flux vector
            I0,Q0,U0,V0 = pmd.read_dap_output(phase,sza,emi,fname)
            
            # append these values to the stokes vectors lists
            # multiply the stokes vectors with the Gaussian weights. 
            Il.append(I0[0]*kdgw[1,i])
            Ql.append(Q0[0]*kdgw[1,i])
            Ul.append(U0[0]*kdgw[1,i])
            Vl.append(V0[0]*kdgw[1,i])
        
        I.append(sum(Il))
        Q.append(sum(Ql))
        U.append(sum(Ul))
        V.append(sum(Vl))
        
    # convert to array
    I = np.array(I)
    Q = np.array(Q)
    U = np.array(U)
    V = np.array(V)
    
    return I,Q,U,V
    
    
def getSPIRflux(lmin=None, lmax=None,idx=0,det=0):
    """
    This function gets the NORMALIZED flux vectors from the SPICAV-IR routine
    within a requested given range of wavelengths.
    You can input an orbital epoch also, default is set at index=0.
    
    Returns: radiance at a given detector,wavelength
    
    Authr: G. Mahapatra
    """
    if det == 0:
        r0 = data.r0[:,idx]
        w0 = data.w0[:,idx]
    else:
        r0 = data.r1[:,idx]
        w0 = data.w1[:,idx]
        
    # drop the "nan's" and normalize the radiance
    w0 = w0[~np.isnan(r0)]
    r0 = r0[~np.isnan(r0)]  # remove nan's
        
    if (lmin==None) & (lmax == None):
        # get the SPICAV-IR radiance and wavelength data
        r = r0
        w = w0
    else:
        r = r0[(w0>lmin)&(w0<lmax)]
        w = w0[(w0>lmin)&(w0<lmax)]
        
    r = r/r[-1]
        
    return r,w
    
def getSPIRpol(lmin=None, lmax=None,idx=0,det=0):  
    """
    This function returns polarized flux data for a given orbit geometry =idx.
    
    It does not do any processing as done by Loic Rossi to arrive at his selected 
    bands.
    
    pol = (det1 - det0)/(det0 + det1)
    
    Returns: pol,wavelength
    
    Authr: G. Mahapatra
    
    """
    
    r0 = data.r0[:,idx]
    w0 = data.w0[:,idx]
    
    r1 = data.r1[:,idx]
    w1 = data.w1[:,idx]
    
    
    if (lmin==None) & (lmax == None):
    # get the SPICAV-IR radiance and wavelength data
        r0 = r0
        w0 = w0
        r1 = r1
        w1 = w1
    else:
        lidx = (w0>lmin)&(w0<lmax)
        r0 = r0[lidx]
        w0 = w0[lidx]
        r1 = r1[lidx]
        w1 = w1[lidx]
        
    # calculate the polarization at all the available geometries
    # polarization is calculated as (det1 - det0)/(det0 + det1)
    pol = np.nan_to_num((r1 - r0)/(r0 + r1))
    
    return pol,w1




#plt.figure(figsize=(12,5))
#style = 'dotted'
#plt.subplot(121);plt.loglog(absTDeep[0:],Iar[0:],c='black',ls=linestyles[style],lw=1)
##plt.subplot(122);plt.semilogx(absTDeep[0:333],Par[0:333],c='black',label='{:3.0f} deg'.format(sza[0]),ls=linestyles[style],lw=1)
#plt.subplot(122);plt.semilogx(absTDeep[0:],Par[0:],c='black',ls=linestyles[style],lw=1,label='With Cloud + Haze at 70 km')
#plt.subplot(121);plt.xlabel('$\sum b_{abs}^{m}$',fontsize='large')
#plt.subplot(122);plt.xlabel('$\sum b_{abs}^{m}$',fontsize='large')
