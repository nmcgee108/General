#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 15:05:56 2026

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
lon_ref, lat_ref = ctd_lons[2],ctd_lats[2] # Innermost station is cast 4 (index 2)
    
ctd_distances = []

for i in range(len(ctd_castnums)):  
    distance = find_distance(ctd_lats[i], ctd_lons[i], lat_ref, lon_ref)
    ctd_distances.append(distance)
    print(distance)


max_dist = max(ctd_distances) 
min_dist = min(ctd_distances) 
dist_colormap = plt.colormaps['viridis'] # Choose colormap
dist_norm = mcolors.Normalize(vmin=min_dist, vmax=max_dist)  # Normalize the colormap with a max and min value
dist_sm = plt.cm.ScalarMappable(cmap=dist_colormap, norm=dist_norm)   # Creates coloring capabilities based on numerical values

# Make plot
fig, axes = plt.subplots(1, 1, figsize=(14, 8))

# Plot data

cast_range = (7, 13)    # Cast 3 (upstream fjord)

nitrate_value = pd.to_numeric(nutrients_file['NO3'][cast_range[0]:cast_range[1]])
sample_cast = nutrients_file['St#'][cast_range[0]:cast_range[1]]
sample_depth = pd.to_numeric(nutrients_file['Depth  '][cast_range[0]:cast_range[1]])

sample_depth_reversed = sample_depth.iloc[::-1]
nitrate_value_reversed = nitrate_value.iloc[::-1]

# Create a linearly interpolated nitrate profile based of cast 8 nitrate data
nitrate_profile = np.interp(ctd_depth, sample_depth_reversed, nitrate_value_reversed)

axes.scatter(nitrate_profile, ctd_temp[1], s = 3)

# for i in range(len(ctd_distances)):  
#     color = dist_colormap(dist_norm(ctd_distances[i]))
#     axes.scatter(ctd_sal[i], ctd_temp[i], color=color, s = 3)

axes.set_xlabel("Nitrate []")
axes.set_ylabel("Conservative Temperature [°C]")
#axes.set_xlim(31.0, 35)
#axes.set_ylim(-2, 3.5)

# axes.legend(loc = "lower right")
# cbar = fig.colorbar(dist_sm, ax=axes, orientation='vertical')
# cbar.set_label("Distance from Inner Station [km]")
















