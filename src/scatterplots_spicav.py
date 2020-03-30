#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 13:42:33 2019

This script creates scatterplots of SPICAV-IR data with gas opacities.
This was a suggestion byy Emmanuel from LATMOS.

@author: gouravmahapatr, TU Delft
"""
import numpy as np
import matplotlib.pyplot as plt
import os
#import linestyles as ls
path = '/Users/gouravmahapatr/Dropbox/PhD/spicav_data'
os.chdir(path)
import read_spicav_ir as rd
from scipy import interpolate

#%% get the total gas opacities as a function of wavelength. 
# we use the absorption opacities derived for the input in line-by-line.
gasop_path = '/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/LineByLine0/'

# load the gas opacities
bmabs_lbl = np.load(gasop_path+'bmabs.npy')

# load the wavelength associated with the gas opacities
wav = np.loadtxt(gasop_path+'wav.dat')

# sum the gas opacities across all the layers
bmabs_lbl_tot = np.sum(bmabs_lbl,axis=1,keepdims=True)

# plot
plt.figure(figsize=[8,6]) 
plt.semilogy(wav,bmabs_lbl_tot,c='k',alpha=0.8)
plt.xlim([1.4,1.5])
plt.xlabel('Wavelength ($\mu m$)',fontsize='x-large')
plt.ylabel('Total gas opacity',fontsize='x-large')
plt.grid()
plt.tight_layout()

#%% get the SPICAV-IR data 
# load the data from a specific orbit
if 'data' not in locals():
    data = rd.processing_orbits(2269,8)

# activate the get_bande()
data.get_bande()

# get the latitude 
lat = data.geo.Lat

# get the spicav wavelength
# Note: The actual wavelength differs b/w two flux channels!!
wav_spicav = data.band_wd0[:,0]

# get the values within 1400 and 1500 nm
idx = np.squeeze(np.array(np.where((wav_spicav<=wav.max()*1.e3)&(wav_spicav>=wav.min()*1.e3))))

wav_spicav_lim = wav_spicav[idx]

#%% get the bmabs_lbl_tot values for the same wavelength points as the SPICAV-IR

# create an interpolation box with bmabs and wavelength information from lbl.
int_bmabs = interpolate.interp1d((wav*1.e3),bmabs_lbl_tot,kind='linear',axis=0)

# get the bmabs_tot corr. to spicav-ir wavelengths
bmabs_spicav = np.squeeze(int_bmabs(wav_spicav_lim))

bmabs_spicav = np.log10(bmabs_spicav)
#%%now load the polarization and flux data
flux = data.band_intensity

pol = data.band_pol

dpol = data.band_dpol

#%% trim them to the wavelengths within limits

flux_lim = flux[idx,:]

pol_lim = pol[idx,:]

dpol_lim = dpol[idx,:]

#%% scatter plots
size = 20
latidx = [0,100,200,300,400,500,-1]

fig = plt.figure(figsize=[15,5])
plt.subplot(121)
ax = plt.gca()
plt.scatter(bmabs_spicav,flux_lim[:,0]/flux_lim[:,0].max(),label='{:2.0f} deg'.format(lat[0]),s=size)
plt.scatter(bmabs_spicav,flux_lim[:,100]/flux_lim[:,100].max(),label='{:2.0f} deg'.format(lat[100]),s=size)
plt.scatter(bmabs_spicav,flux_lim[:,200]/flux_lim[:,200].max(),label='{:2.0f} deg'.format(lat[200]),s=size)
plt.scatter(bmabs_spicav,flux_lim[:,300]/flux_lim[:,300].max(),label='{:2.0f} deg'.format(lat[300]),s=size)
plt.scatter(bmabs_spicav,flux_lim[:,400]/flux_lim[:,400].max(),label='{:2.0f} deg'.format(lat[400]),s=size)
plt.scatter(bmabs_spicav,flux_lim[:,500]/flux_lim[:,500].max(),label='{:2.0f} deg'.format(lat[500]),s=size)
plt.scatter(bmabs_spicav,flux_lim[:,-3]/flux_lim[:,-3].max(),label='{:2.0f} deg'.format(lat[-3]),s=size)
plt.legend()
plt.xlabel('log(Total gas opacity)',fontsize='x-large')
plt.ylabel('Total flux',fontsize='x-large')
plt.grid()
#plt.xlim([0,300])

plt.subplot(122)
plt.scatter(bmabs_spicav,pol_lim[:,0],label='{:2.0f} deg'.format(lat[0]),s=size)
plt.scatter(bmabs_spicav,pol_lim[:,100],label='{:2.0f} deg'.format(lat[100]),s=size)
plt.scatter(bmabs_spicav,pol_lim[:,200],label='{:2.0f} deg'.format(lat[200]),s=size)
plt.scatter(bmabs_spicav,pol_lim[:,300],label='{:2.0f} deg'.format(lat[300]),s=size)
plt.scatter(bmabs_spicav,pol_lim[:,400],label='{:2.0f} deg'.format(lat[400]),s=size)
plt.scatter(bmabs_spicav,pol_lim[:,500],label='{:2.0f} deg'.format(lat[500]),s=size)
plt.scatter(bmabs_spicav,pol_lim[:,-3],label='{:2.0f} deg'.format(lat[-3]),s=size)
plt.legend()
plt.xlabel('log(Total gas opacity)',fontsize='x-large')
plt.ylabel('Polarization',fontsize='x-large')
plt.grid()
#plt.xlim([0,300])
plt.tight_layout()

#%% create scatter plots with polarization of model data and observation data
# load the wavelength
wav_model = np.load('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/wav.npy')

# load the points of interest
wavscat = np.array([1,65,159])

# load the bmabs corresponding to these values
bmabs_model = np.log10(int_bmabs(wav_model[wavscat]*1.e3))

#%% load the model data
j = 17

# 
kspectra_ct65 = np.load('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/geos_source/geosOutcloud+hazeext/kspectra_ct65.npy')
kspectra_ct70 = np.load('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/geos_source/geosOutcloud+hazeext/kspectra_ct70.npy')
kspectra_ct75 = np.load('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/fullKDIS/geos_source/geosOutcloud+hazeext/kspectra_ct75.npy')


Ic65 = kspectra_ct65[j,0,wavscat]
Qc65 = kspectra_ct65[j,1,wavscat]
Uc65 = kspectra_ct65[j,2,wavscat]
Pc65 = -Qc65/Ic65

Ic70 = kspectra_ct70[j,0,wavscat]
Qc70 = kspectra_ct70[j,1,wavscat]
Uc70 = kspectra_ct70[j,2,wavscat]
Pc70 = -Qc70/Ic70

Ic75 = kspectra_ct75[j,0,wavscat]
Qc75 = kspectra_ct75[j,1,wavscat]
Uc75 = kspectra_ct75[j,2,wavscat]
Pc75 = -Qc75/Ic75

# %% plot the model and observation data
#fig,ax = plt.subplots(1,1,figsize=[10,7])
#colors = ['blue','black','green','blue','green']

plt.figure(figsize=[15,5])
ax1 = plt.subplot(122)



size = 20
#ax.scatter(bmabs_spicav,pol_lim[:,100],label='{:2.0f} deg'.format(lat[100]),s=size)
#ax.scatter(bmabs_spicav,pol_lim[:,200],label='{:2.0f} deg'.format(lat[200]),s=size)
#ax.scatter(bmabs_spicav,pol_lim[:,300],label='{:2.0f} deg'.format(lat[300]),s=size)
#ax.scatter(bmabs_spicav,pol_lim[:,400],label='{:2.0f} deg'.format(lat[400]),s=size)
ax1.scatter(bmabs_spicav,pol_lim[:,500],label='(data) lat = {:2.0f} deg'.format(lat[500]),s=size,alpha=0.3)
ax1.scatter(bmabs_spicav,pol_lim[:,-3],label='(data) lat = {:2.0f} deg'.format(lat[-3]),s=size,alpha=0.3)
ax1.legend()

#%%
size = 150
for i in range(len(wavscat)):
    if i == 1: 
        ax1.scatter(bmabs_model[i],Pc65[i],marker='o',c='g',s=size,label='(model) Cloud at 65 km, lat = {:2.0f}'.format(j))
        ax1.scatter(bmabs_model[i],Pc70[i],marker='*',c='g',s=size,label='(model) Cloud at 70 km, lat = {:2.0f}'.format(j))
        ax1.scatter(bmabs_model[i],Pc75[i],marker='v',c='g',s=size,label='(model) Cloud at 75 km, lat = {:2.0f}'.format(j))
    else:
        ax1.scatter(bmabs_model[i],Pc65[i],marker='o',c='g',s=size)
        ax1.scatter(bmabs_model[i],Pc70[i],marker='*',c='g',s=size)
        ax1.scatter(bmabs_model[i],Pc75[i],marker='v',c='g',s=size) 
ax1.legend(loc=3,fontsize='medium')
ax1.set_xlabel('log(Total gas opacity)', fontsize='x-large')
ax1.set_ylabel('Polarization', fontsize='x-large')

#%% flux plots

ax2 = plt.subplot(121)


size = 20
#ax.scatter(bmabs_spicav,pol_lim[:,100],label='{:2.0f} deg'.format(lat[100]),s=size)
#ax.scatter(bmabs_spicav,pol_lim[:,200],label='{:2.0f} deg'.format(lat[200]),s=size)
#ax.scatter(bmabs_spicav,pol_lim[:,300],label='{:2.0f} deg'.format(lat[300]),s=size)
#ax.scatter(bmabs_spicav,pol_lim[:,400],label='{:2.0f} deg'.format(lat[400]),s=size)
ax2.scatter(bmabs_spicav,flux_lim[:,500]/flux_lim[:,500].max(),label='(data) lat = {:2.0f} deg'.format(lat[500]),s=size,alpha=0.3)
ax2.scatter(bmabs_spicav,flux_lim[:,-3]/flux_lim[:,-3].max(),label='(data) lat = {:2.0f} deg'.format(lat[-3]),s=size,alpha=0.3)
ax2.legend()

#%%
size = 150
for i in range(len(wavscat)):
    if i == 1: 
        ax2.scatter(bmabs_model[i],Ic65[i]/Ic65[0].max(),marker='o',c='g',s=size,label='(model) Cloud + haze at 65 km, lat = {:2.0f}'.format(j))
        ax2.scatter(bmabs_model[i],Ic70[i]/Ic70[0].max(),marker='*',c='g',s=size,label='(model) Cloud + haze at 70 km, lat = {:2.0f}'.format(j))
        ax2.scatter(bmabs_model[i],Ic75[i]/Ic75[0].max(),marker='v',c='g',s=size,label='(model) Cloud + haze at 75 km, lat = {:2.0f}'.format(j))
    else:
        ax2.scatter(bmabs_model[i],Ic65[i]/Ic65[0].max(),marker='o',c='g',s=size)
        ax2.scatter(bmabs_model[i],Ic70[i]/Ic70[0].max(),marker='*',c='g',s=size)
        ax2.scatter(bmabs_model[i],Ic75[i]/Ic75[0].max(),marker='v',c='g',s=size) 
ax2.legend(loc=3,fontsize='medium')
ax2.set_xlabel('log(Total gas opacity)', fontsize='x-large')
ax2.set_ylabel('Flux', fontsize='x-large')