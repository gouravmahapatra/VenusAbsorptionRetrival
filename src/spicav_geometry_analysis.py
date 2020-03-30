#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 11:16:33 2018

This script reads the SPICAV IR data from the "read_spicav_ir.py" script
made by Loic Rossi and does some observation geometry data plotting.

@author: gouravmahapatr, TU Delft
"""
import os
import numpy as np
import matplotlib.pyplot as plt
# first load all the geometry information related to a particular orbit from the VEX data
path = '/Users/gouravmahapatr/Dropbox/PhD/spicav_data'
os.chdir(path)
import read_spicav_ir as rd
import pymiedap.pymiedap as pmd

# get the orbits list that are within the orbit range
orbit_list = rd.find_orbits(400,1500)


    
#data = rd.processing_orbits(100,2)
#localtime = data.geo.local_time
#lat = data.geo.Lat
#phase = data.geo.phase
#sza = data.geo.sza
#emission = data.geo.emission
localtime = []
lat = []
phase = []
sza = []
emi = []
azi = []

for i in range(len(orbit_list)):
    orbit_n = orbit_list[i,0]
    orbit_a = orbit_list[i,1]
    data = rd.processing_orbits(orbit_n,orbit_a)
    data.geo.calc_azimuth()
    localtime.extend(list(data.geo.local_time))
    lat.extend(list(data.geo.Lat))
    phase.extend(list(data.geo.phase))
    sza.extend(list(data.geo.sza))
    emi.extend(list(data.geo.emission))
    azi.extend(list(data.geo.azimuth))
    

cosbeta = pmd.get_cosbeta(np.array(phase),np.array(sza),np.array(emi),np.array(azi)) # in radians

beta = np.degrees(np.arccos(cosbeta)) # Converted into degrees


 

''' Make the geometries file that will be used as an input to the geos_code'''

# first make sure the angles are all in array
alpha = np.array(phase)
sza = np.array(sza)
emi = np.array(emi)
azi = np.array(azi)
beta = np.array(beta)
localtime = np.array(localtime)
lat = np.array(lat)

# start making data manipulations

# only keep those values that have lat range of [-10,10] deg
alpha1 = alpha[(lat<5)&(lat>-5)]
sza1 = sza[(lat<5)&(lat>-5)]
emi1 = emi[(lat<5)&(lat>-5)]
azi1 = azi[(lat<5)&(lat>-5)]
beta1 = beta[(lat<5)&(lat>-5)]
localtime1 = localtime[(lat<5)&(lat>-5)]
lat1 = lat[(lat<5)&(lat>-5)]

# now remove all the values out of the daylight zone
alpha1 = alpha1[(localtime1>6)&(localtime1<18)]
sza1 = sza1[(localtime1>6)&(localtime1<18)]
emi1 = emi1[(localtime1>6)&(localtime1<18)]
azi1 = azi1[(localtime1>6)&(localtime1<18)]
beta1 = beta1[(localtime1>6)&(localtime1<18)]
lat1 = lat1[(localtime1>6)&(localtime1<18)]
localtime1 = localtime1[(localtime1>6)&(localtime1<18)]

# rearrange them as per local time bins
# only take one value out of the whole sorted bin values
lt = localtime1 
alpha2 =[]
sza2 =[]
emi2 =[]
azi2 = []
beta2 = []
lat2 = []
localtime2 = []
nvals,edges = np.histogram(localtime1,bins=100)
for i in range(len(edges)):
    if alpha1[(localtime1>edges[i])&(localtime1<[i+1])].size > 0: 
        alpha2.append(alpha1[(localtime1>edges[i])&(localtime1<[i+1])][0])
        sza2.append(sza1[(localtime1>edges[i])&(localtime1<[i+1])][0])
        emi2.append(emi1[(localtime1>edges[i])&(localtime1<[i+1])][0])
        azi2.append(azi1[(localtime1>edges[i])&(localtime1<[i+1])][0])
        beta2.append(beta1[(localtime1>edges[i])&(localtime1<[i+1])][0])
        lat2.append(lat1[(localtime1>edges[i])&(localtime1<[i+1])][0])
        localtime2.append(localtime1[(localtime1>edges[i])&(localtime1<[i+1])][0])

# set all the beta's to zero for now
beta2 = np.zeros(beta.shape)

# perform the writing
ngeos = len(alpha2)
geosfile = open('/Users/gouravmahapatr/Dropbox/PhD/Codes/VenusAbsCode/geos.in','w')
geosfile.write('# this file contains the geometries obtained from VEx SPICAV IR dataset \n')
geosfile.write('# number of available geometries (ngeos): \n')
geosfile.write(str(ngeos)+'\n')
geosfile.write('#    alpha    theta0    theta    phi    beta \n')
for i in range(ngeos):
    geosfile.write('    '+str('{:5.3f}').format(alpha2[i])+'    ')
    geosfile.write(str('{:6.2f}').format(sza2[i])+'    ')
    geosfile.write(str('{:6.2f}').format(emi2[i])+'    ')
    geosfile.write(str('{:6.2f}').format(azi2[i])+'    ')
    geosfile.write(str('{:6.2f}').format(beta2[i])+'    \n')
geosfile.close()
"""
-------------------------------------------------------------------------------
"""
# make plots
#os.chdir('../figures/')
#plt.scatter(localtime,lat,c=emi)
#plt.xlim([6,18])
#plt.ylim([-80,80])
#clb = plt.colorbar()
#clb.set_label('Emission angle (deg)',fontsize='large')
#plt.xlabel('Local time',fontsize='large')
#plt.ylabel('Latitude (deg)',fontsize='large')
#plt.tight_layout()
#plt.savefig('emission_plot_400to1500.pdf')
#plt.close()

#plt.scatter(localtime,lat,c=sza)
#plt.xlim([6,18])
#plt.ylim([-80,80])
#clb = plt.colorbar()
#clb.set_label('Solar zenith angle (deg)',fontsize='large')
#plt.xlabel('Local time',fontsize='large')
#plt.ylabel('Latitude (deg)',fontsize='large')
#plt.tight_layout()
#plt.savefig('sza_plot_400to1500.pdf')
#plt.close()
#
#plt.scatter(localtime,lat,c=beta)
#plt.xlim([6,18])
#plt.ylim([-80,80])
#clb = plt.colorbar()
#clb.set_label('Rotation angle (deg)',fontsize='large')
#plt.xlabel('Local time',fontsize='large')
#plt.ylabel('Latitude (deg)',fontsize='large')
#plt.tight_layout()
#plt.savefig('beta_plot_400to1500.pdf')
#plt.close()

