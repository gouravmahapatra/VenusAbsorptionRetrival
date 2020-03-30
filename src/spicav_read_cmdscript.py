#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 11:32:09 2018

Script to utilize the read_spicav_ir.py function to read and plot spicav IR data set.

@author: Loic Rossi, gouravmahapatr
"""

import read_spicav_ir as rd
import matplotlib.pyplot as plt

"""
load a specific orbit on the basis of its orbit number 
"""
data = rd.processing_orbits(1478,8)


plt.plot(data.geo.Alt_sc)  # plt the spacecraft altitude 
plt.plot(data.geo.time)    # plot the spacecraft time stamp associated with geometry

data.geo.time.shape

plt.plot(data.geo.phase, data.geo.Alt_sc) # plot the phase angle vs. SC alt
plt.plot(data.geo.phase, data.geo.Lat)    # plot phase angle vs. latitude

plt.plot(data.geo.phase, data.geo.emission)
plt.plot(data.geo.phase, data.geo.sza)
data.w0
data.w0.shape
data.w1.shape
data.w0
data.w1
data.r0
data.r1
data.r1.shape
data.wd0
data.wd0.shape
data?
data.num_orbit
data.num_orbit.shape
plt.scatter(data.wd0[:,0],data.rd0[:,0])
plt.scatter(data.wd0[:,10],data.rd0[:,10])
plt.scatter(data.wd1[:,10],data.rd1[:,10])
plt.scatter(data.w1[:,10],data.r1[:,10])
plt.scatter(data.phase, data.pol[0,:])
plt.scatter(data.phase, data.pol[4,:])
data.wd0[:,0]
plt.scatter(data.phase, data.pol[8,:])
plt.scatter(data.geo.Alt_sc, data.pol[8,:])
data
data = rd.processing_orbits(1000,0,orbit_max=1200,orbit_range=True)
data = rd.processing_orbits([1478,2048],[8,8])
data = rd.processing_orbits([1478,2048],[8,8],list=True)
data = rd.processing_orbits([1478,2269],[8,8],list=True)
data
dataA = rd.processing_orbits(1478,8)
dataB = rd.processing_orbits(2269,8)
data = dataA + dataB
data
dataF = data.filter(criterium='latitude',vmin=0,vmax=30)
dataF
data.phase.min()
dataF.phase.min()
dataF.geo.Lat.min()
dataF.geo.Lat.max()
data.geo.Lat.max()
data.filter?
data.get_bande?
data.get_bande()
data.show_band?
data.show_band()

data = rd.processing_orbits([1478,2269],[8,8],list=True)
data = rd.processing_orbits(1478,8)
.ge
data.get_bande()
plot(data.phase[34,:],data.band_pol[34,:])
plt.scatter(data.phase[34,:],data.band_pol[34,:])
plt.scatter(data.phase,data.band_pol[34,:])
plt.scatter(data.band_wd0[:,20],data.band_pol[:,20])
plt.scatter(data.band_wd0[:,20],data.band_pol[:,20],'+')
plt.scatter(data.band_wd0[:,20],data.band_pol[:,20],ms=3)
plt.scatter(data.band_wd0[:,20],data.band_pol[:,20],s=3)
plt.scatter(data.band_wd0[:,20],data.band_pol[:,45],s=3)
plt.scatter(data.band_wd0[:,20],data.band_pol[:,100],s=3)
data.phase
data.phase[-10]
data.phase[-15]
plt.scatter(data.band_wd0[:,-15],data.band_pol[:,-15],s=3)
data.phase_pol_colat(4)
data.show_band()