#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 18:48:35 2019

This script creates errorbars for the degree of polarization 
of SPICAV-IR orbits.

@author: gouravmahapatr, TU Delft
"""

import matplotlib.pyplot as plt
import sys
import os
path = '/Users/gouravmahapatr/Dropbox/PhD/spicav_data'
os.chdir(path)
import read_spicav_ir as rd
import linestyles as ls


# load the data from a specific orbit
data = rd.processing_orbits(1478,8)


# activate the get_bande()
data.get_bande()


# define the data step
step = 20

# define the index
idx = 0

# get the latitude 
lat = data.geo.Lat[idx]

pol = data.band_pol[::step,idx] 
pol+=0.0

# plot the errorbars 
plt.errorbar(data.band_wd0[::step,idx],pol,yerr=data.band_dpol[::step,idx],label=lat,c='k',ls=ls.linestyles['densely dashdotted'])

plt.legend(title='Latitude (deg)')