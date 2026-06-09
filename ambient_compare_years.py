#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 16:30:43 2026

@author: nataliemcgee
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
import gsw



# Load the CTD NetCDF dataset
ctd_file = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Final CTD Datasets/uc_patch_dataset_new.nc"

with Dataset(ctd_file, 'r') as ncfile:
    
    # Extract cast and depth-indexed variables (many values per cast)
    # Creates a 2D numpy array of values indexed by cast number and depth bin (missing data filled with NaN)
    depth_bin = np.ma.filled(ncfile.variables['Depth Bin'][:], np.nan)
    SA = np.ma.filled(ncfile.variables['Absolute Salinity'][:], np.nan)
    CT = np.ma.filled(ncfile.variables['Conservative Temp'][:], np.nan)

file_list = [ "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Morven CTD Data 2013-2019/Morven Data - 2013.csv",
             "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Morven CTD Data 2013-2019/Morven Data - 2015.csv",
             "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Morven CTD Data 2013-2019/Morven Data - 2016.csv",
             "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Morven CTD Data 2013-2019/Morven Data - 2019.csv"]


colors = plt.cm.managua(np.linspace(0, 1, 5))
years = ["2013", "2015", "2016", "2019", "2024"]

#-----------------------------
# THINGS TO ENTER
#-----------------------------
cast_num = 3
    
lat = 72.9167
zi = -depth_bin[cast_num-1][::-1] # Must input values from deepest to shallowest!
Ta_2024_rough = CT[cast_num-1][::-1]   # ambient fjord temperature at zi
Sa_2024_rough = SA[cast_num-1][::-1] # ambient fjord salinity at zi
p_2024= gsw.p_from_z(zi,lat)


Ta_2024= savgol_filter(Ta_2024_rough, window_length=10, polyorder=3)
Sa_2024= savgol_filter(Sa_2024_rough, window_length=10, polyorder=3)

n_squared_2024_rough,_ = gsw.Nsquared(Sa_2024, Ta_2024, p_2024, lat)

n_squared_2024= savgol_filter(n_squared_2024_rough[50:], window_length=80, polyorder=1)

fig, axes = plt.subplots(1, 3, figsize=(11, 7), sharey=True)

for i in range(len(file_list)):
    
    df = pd.read_csv(file_list[i])
    z = -df["z"]
    Ta_rough = df["Ta"]
    Sa_rough = df["Sa"]
    p = gsw.p_from_z(z,lat)
    
    Ta= savgol_filter(Ta_rough, window_length=10, polyorder=3)
    Sa= savgol_filter(Sa_rough, window_length=10, polyorder=3)
    
    n_squared_rough,_ = gsw.Nsquared(Sa, Ta, p, lat)
    n_squared= savgol_filter(n_squared_rough[50:], window_length=80, polyorder=1)
    
    
    axes[0].plot(Ta, z, color = colors[i], label = years[i])
    axes[1].plot(Sa, z, color = colors[i], label = years[i])
    axes[2].plot(n_squared, z[51:], color = colors[i], label = years[i]) # N_squared has one fewer entry than the others
                                                                         # thus it is offset by 0.5m to do it this way but thats ok
   
axes[0].plot(Ta_2024, zi, color = "darkblue", linewidth = 2.75, label = years[-1])
axes[1].plot(Sa_2024, zi, color = "darkblue", linewidth = 2.75, label = years[-1])
axes[2].plot(n_squared_2024, zi[51:], color = "darkblue", linewidth = 2.75, label = years[-1])

#axes[2].set_xscale("log")

axes[0].set_ylabel("Depth [m]", fontsize = 13)
axes[0].set_xlabel("CT [°C]", fontsize = 13)
axes[1].set_xlabel(r"S$_A$ [g/kg]", fontsize = 13)
axes[2].set_xlabel(r"N$^2$ [s$^{-2}$]", fontsize = 13)

axes[0].set_ylim(-750, 0)
axes[0].set_xlim(0, 3.5)
axes[1].set_xlim(33.5, 35.1)
axes[2].set_xlim(-5e-6,3.5e-5)

axes[0].grid(True)
axes[1].grid(True)
axes[2].grid(True)
    
axes[1].legend(fontsize = 14)




plt.suptitle("Ambient Temperature and Salinity Stratifications")
plt.tight_layout()
plt.show()