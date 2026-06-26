#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 11:31:49 2026

@author: nataliemcgee
"""

import xarray as xr
import numpy as np
import gsw
import matplotlib.pyplot as plt
import math
import pandas as pd
import matplotlib.colors as mcolors

plt.rcParams['font.size']=14

ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"
nutrients_file = pd.read_csv("/Users/nataliemcgee/Documents/Upernavik Data/Nutrients/NutrientsUS2024_plotting.csv", encoding="latin-1")

ctd_ds = xr.open_dataset(ctd_netcdf)

cast_range = (2, 9)  # Enter casts of interest

ctd_depth = ctd_ds["depth"].values
ctd_sal = ctd_ds["SAL_ABSOLUTE"][cast_range[0]-1:cast_range[1]].values
ctd_temp = ctd_ds["CONSERVATIVE_TEMP"][cast_range[0]-1:cast_range[1]].values
ctd_lats = ctd_ds["LAT"][cast_range[0]-1:cast_range[1]].values
ctd_lons = ctd_ds["LON"][cast_range[0]-1:cast_range[1]].values
ctd_castnums = ctd_ds["cast"][cast_range[0]-1:cast_range[1]].values +1

sigma0 = gsw.sigma0(ctd_sal, ctd_temp)

# Extract nutrient data
nitrate_value = pd.to_numeric(nutrients_file['NO3'][1:41])
sample_cast = nutrients_file['St#'][1:41]
sample_depth = -pd.to_numeric(nutrients_file['Depth  '][1:41])

def find_distance(lat1, lon1, lat2, lon2):

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

# Calculate the distance of each cast from the innermost station
lon_ref, lat_ref = ctd_lons[3],ctd_lats[3] # Innermost station is cast 4 (index 3)
    
ctd_distances = []

for i in range(len(ctd_castnums)):   
    ctd_distances.append(find_distance(ctd_lats[i], ctd_lons[i], lat_ref, lon_ref))


max_dist = max(ctd_distances) 
min_dist = min(ctd_distances) 
dist_colormap = plt.colormaps['viridis'] # Choose colormap
dist_norm = mcolors.Normalize(vmin=min_dist, vmax=max_dist)  # Normalize the colormap with a max and min value
dist_sm = plt.cm.ScalarMappable(cmap=dist_colormap, norm=dist_norm)   # Creates coloring capabilities based on numerical values

# Make plot
fig, axes = plt.subplots(1, 1, figsize=(14, 8))

# Create grid of salinity and temperature
minS = 31
maxS = 35
minT = -2
maxT = 4

sx = np.arange(minS, maxS + 0.1, 0.1)
ty = np.arange(minT, maxT + 0.1, 0.1)

S, T = np.meshgrid(sx, ty)
PDEN = gsw.rho(S, T, 0) - 1000  # Potential density anomaly at reference pressure = 0

# Manually label contours
contour = axes.contour(S, T, PDEN, levels=[24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28], colors='grey')
label_positions = [(31.4, -1), (31.9, -1), (32.4, -1), (33.2, 3.1), (33.8, 3.1), (34.5, 1.4), (35.0, 1)]
plt.clabel(contour, inline=True, manual=label_positions, fmt='%1.1f')

# Mixing line params
T_deepwater = 2.550228
S_deepwater = 34.697933
T_ice_eff = -90
    
def runoff_line(x):
    T_deepwater = 2.550228
    S_deepwater = 34.697933
    slope = T_deepwater/S_deepwater
    return slope*x

def melting_line(x):
    slope = (T_deepwater-T_ice_eff)/S_deepwater
    return slope*x+T_ice_eff

axes.plot(sx, runoff_line(sx), color ="orange", linestyle = "dashed", 
             label = "Runoff Line", 
             zorder = 0)
axes.plot(sx, melting_line(sx), color ="pink", linestyle = "dashed", 
             label = "Melting Line", 
             zorder = 0)

# Plot TS data

for i in range(len(ctd_distances)):  
    color = dist_colormap(dist_norm(ctd_distances[i]))
    axes.scatter(ctd_sal[i], ctd_temp[i], color=color, s = 3)

axes.set_xlabel("Absolute Salinity [g/kg]")
axes.set_ylabel("Conservative Temperature [°C]")
axes.set_xlim(31.0, 35)
axes.set_ylim(-2, 3.5)

axes.legend(loc = "lower right")
cbar = fig.colorbar(dist_sm, ax=axes, orientation='vertical')
cbar.set_label("Distance from Inner Station [km]")

# Plot nutrient data

colormap = plt.colormaps['YlGn']
norm = mcolors.Normalize(vmin=0, vmax=18)
sm_nitrate = plt.cm.ScalarMappable(cmap=colormap, norm=norm)

for j in range(len(nitrate_value), 0, -1): #plot shallow samples first
    index = np.where(ctd_castnums==sample_cast[j])[0]
    index = int(index[0])
    sal = ctd_sal[index][-sample_depth[j]]
    temp = ctd_temp[index][-sample_depth[j]]
    #color = colormap(norm(nitrate_value[j]))
    #axes.scatter(sal, temp, color=color, edgecolor = 'k', s= 80, zorder=4)
    

    color = dist_colormap(dist_norm(ctd_distances[index]))
    axes.scatter(sal, temp, color=color, marker = "o", edgecolor = 'k', s=nitrate_value[j]*15, zorder=4)
    
# cbar_ax = fig.add_axes([0.87,0.1,0.021,0.78])
# cbar = fig.colorbar(sm_nitrate, cax=cbar_ax)
# cbar.set_label('Nitrate [µmol N-NO3/L]')






















