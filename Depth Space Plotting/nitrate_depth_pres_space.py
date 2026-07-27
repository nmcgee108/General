#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 11:09:44 2026

@author: nataliemcgee
"""

import xarray as xr
import numpy as np
import gsw
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from matplotlib.colors import PowerNorm
import math
import pandas as pd
import matplotlib.colors as mcolors

plt.rcParams['font.size']=18

nutrient_profiles = "/Users/nataliemcgee/Documents/Upernavik Data/Nutrients/nitrate_profiles.nc"
ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"


nitrate_ds = xr.open_dataset(nutrient_profiles)
ctd_ds = xr.open_dataset(ctd_netcdf)

ctd_sal = ctd_ds["SAL_ABSOLUTE"].values
ctd_temp = ctd_ds["CONSERVATIVE_TEMP"].values
ctd_castnums = ctd_ds["cast"].values

sigma0 = gsw.sigma0(ctd_sal, ctd_temp)


depth = nitrate_ds["depth"].values
nitrate = nitrate_ds["Nitrate"].values
nitr_lats = nitrate_ds["latitude"].values
nitr_lons = nitrate_ds["longitude"].values
castnums = nitrate_ds["cast"].values

# Pad the ctd data so all casts are the same length as nitrate cast data
target_len = nitrate.shape[1]

sal_padded = np.pad(ctd_sal, ((0, 0), (0, target_len - ctd_sal.shape[1])), constant_values=np.nan)
temp_padded = np.pad(ctd_temp, ((0, 0), (0, target_len - ctd_sal.shape[1])), constant_values=np.nan)
sigma0_padded = np.pad(sigma0, ((0, 0), (0, target_len - sigma0.shape[1])), constant_values=np.nan)


fig, axes = plt.subplots(1, 1, figsize=(8, 10))

def find_distance(lat2, lon2):
    '''finds distance in km between a point and station 4 which is closest to the glacier using
    the Haversine formula'''
    
    # lat and lon of station 4:
    lat1 = nitr_lats[3]
    lon1 = nitr_lons[3]
        
    R = 6371.0  # Earth's radius in km

    # Convert degrees to radians
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    # Haversine formula
    a = math.sin(dphi / 2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    distance = R * c
    return distance

max_dist = find_distance(nitr_lats[7], nitr_lons[7]) #using cast 8

colormap = plt.colormaps['viridis']
norm = mcolors.Normalize(vmin=0, vmax=max_dist)
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)

for i in range(19): #[1, 2, 3, 7, 8, 12, 13, 14, 15, 16, 17, 18]:
    if np.isin(castnums[i], [1, 4, 5, 10, 11, 12]):
        continue
    
    zorder = 1
    linewidth = 2
    label = None
    
    if castnums[i]<9: 
        dist = find_distance(nitr_lats[i], nitr_lons[i])
        color = colormap(norm(dist))
        zorder = 2
        if castnums[i]==3: label = "Fjord casts 2 and 3"
        
    if castnums[i]==9: 
        color = "k"
        zorder = 2
        linewidth = 3
        label = "Cast 9"
        
    if castnums[i]>9: 
        color = "gray"
        linewidth = 1
        if castnums[i]==14: label = "Trough"
    
    axes.plot(nitrate[i], -depth, color = color, zorder = zorder, linewidth=linewidth, label = label)
    #axes.plot(nitrate[i], sigma0_padded[i], color = color, zorder = zorder, linewidth =linewidth, label = label)

axes.legend(loc="lower left")
axes.set_xlabel(r"Nitrate [$\mu$M]")

# DEPTH PROFILE STUFF
axes.set_ylabel(r"Depth [m]")
axes.set_ylim(-500, 50)

# DENSITY PROFILE STUFF
# axes.invert_yaxis()
# axes.set_ylabel(r"Potential Density Anomaly [kg/m$^3$]")
# axes.set_ylim(27.8, 26.25)
# axes.set_xlim(10, 21)
# axes.axhspan(27.05, 27.45, color = 'skyblue', alpha = 0.3)
# axes.axhspan(25.7, 26.8, color = 'thistle', alpha = 0.3)



cbar_ax = fig.add_axes([1.03,0.1,0.05,0.8])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Fjord Profiles: Distance from Innermost Station [km]', fontsize = 14)






