#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 12:40:53 2018

This script reads spicav data and produces some useful output 
for Venus cloud altimetry.

@author: gouravmahapatr, TU Delft
"""
import numpy as np
import os
import matplotlib.pyplot as plt

path = '/Users/gouravmahapatr/Dropbox/PhD/spicav_data'
os.chdir(path)
import read_spicav_ir as rd

"""
load a specific orbit on the basis of its orbit number 
"""
data = rd.processing_orbits(1478,8)

# read the orbit geometries for the pass
phase = data.geo.phase
emission = data.geo.emission
sza = data.geo.sza

plt.plot(phase,sza)