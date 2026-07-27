#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 15:40:17 2026

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

ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Morven CTD Data 2013-2019/2015_profiles.nc"

ctd_ds = xr.open_dataset(ctd_netcdf)

ctd_depth = ctd_ds["depth"].values
ctd_sal = ctd_ds["Absolute_Salinity"].values
ctd_temp = ctd_ds["Conservative_Temperature"].values
ctd_lats = ctd_ds["latitude"].values
ctd_lons = ctd_ds["longitude"].values
ctd_castnums = np.arange(len(ctd_lats))

df = pd.DataFrame({
    "cast": ctd_castnums,
    "latitude": ctd_lats,
    "longitude": ctd_lons
})

df.to_csv("2015_locations.csv", index=False)
print("locations saved")

sigma0 = gsw.sigma0(ctd_sal, ctd_temp)

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
lon_ref, lat_ref = ctd_lons[12],ctd_lats[12] # Innermost station is ?
    
ctd_distances = []

for i in range(len(ctd_castnums)):  
    distance = find_distance(ctd_lats[i], ctd_lons[i], lat_ref, lon_ref)
    ctd_distances.append(distance)

max_dist = max(ctd_distances) 
max_dist = 47
min_dist = min(ctd_distances) 
dist_colormap = plt.colormaps['viridis'] # Choose colormap
dist_norm = mcolors.Normalize(vmin=min_dist, vmax=max_dist)  # Normalize the colormap with a max and min value
dist_sm = plt.cm.ScalarMappable(cmap=dist_colormap, norm=dist_norm)   # Creates coloring capabilities based on numerical values

# Make plot
fig, axes = plt.subplots(1, 1, figsize=(12, 10))

# Create grid of salinity and temperature
minS = 30
maxS = 35
minT = -2
maxT = 8

sx = np.arange(minS, maxS, 0.1)
ty = np.arange(minT, maxT, 0.1)

S, T = np.meshgrid(sx, ty)
PDEN = gsw.rho(S, T, 0) - 1000  # Potential density anomaly at reference pressure = 0

# Manually label contours
contour = axes.contour(S, T, PDEN, levels=[23, 23.5, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28], colors='grey')
label_positions = [(30.5, 4), (30.5, 2), (31.2, 1.5), (32, 2), (32.6, 2), (33.2, 2.1), (33.8, 2.2), (34.5, 1), (35.0, 2)]
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

# FYI: the cast number is the same as the index
for i in [7, 8, 9, 11]:  
    color = dist_colormap(dist_norm(ctd_distances[i]))
    color = "C1"
    if i ==7:
        label = "S branch (7, 8, 9, 11)"
    else: label = ""
    #axes.scatter(ctd_sal[i], ctd_temp[i], color=color, s = 5, label = label)
    
    
for i in [6, 12, 21, 5, 4, 22, 23, 3]:  
    color = dist_colormap(dist_norm(ctd_distances[i]))
    #color = "C0"
    if i ==6:
        label = ""
    else: label = ""
    #axes.scatter(ctd_sal[i], ctd_temp[i], color=color, s = 5, label = label)
    

for i in [16, 17]:  
    color = dist_colormap(dist_norm(ctd_distances[i]))
    color = "C2"
    if i ==16:
        label = "N branch (16, 17)"
    else: label = ""
    #axes.scatter(ctd_sal[i], ctd_temp[i], color=color, s = 5, label = label)
    
for i in [13, 14, 15, 18, 20]:  
    color = dist_colormap(dist_norm(ctd_distances[i]))
    color = "indigo"
    if i ==13:
        label = "Shallows (13, 14, 15, 18, 20)"
    else: label = ""
    axes.scatter(ctd_sal[i], ctd_temp[i], color=color, s = 5, label = label)
    
for i in [0, 1, 2, 29]:  
    color = dist_colormap(dist_norm(ctd_distances[i]))
    color = "C3"
    if i ==0:
        label = "Shelf (0, 1, 2, 29)"
    else: label = ""
    #axes.scatter(ctd_sal[i], ctd_temp[i], color=color, s = 5, label = label)
    
for i in [24, 25, 26, 27, 28]:  
    color = dist_colormap(dist_norm(ctd_distances[i]))
    color = "C4"
    if i ==24:
        label = "NW branch (24, 25, 26, 27, 28)"
    else: label = ""
    axes.scatter(ctd_sal[i], ctd_temp[i], color=color, s = 5, label = label)
    
    
    
axes.set_xlabel("Absolute Salinity [g/kg]")
axes.set_ylabel("Conservative Temperature [°C]")
axes.set_xlim(32, 35)
axes.set_ylim(-1, 3.3)


axes.legend()
cbar = fig.colorbar(dist_sm, ax=axes, orientation='vertical')
cbar.set_label("Distance from Station 12 [km]")






















