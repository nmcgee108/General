#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 10:58:11 2026

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

plt.rcParams['font.size']=14

ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"
bathymetry_file = pd.read_csv("/Users/nataliemcgee/Documents/Upernavik Data/Bathymetry Data/trough_section_bathymetry1.csv")
nutrients_file = pd.read_csv("/Users/nataliemcgee/Documents/Upernavik Data/Nutrients/NutrientsUS2024_plotting.csv", encoding="latin-1")
nutrient_profiles = "/Users/nataliemcgee/Documents/Upernavik Data/Nutrients/nitrate_profiles.nc"

nitrate_ds = xr.open_dataset(nutrient_profiles)
ctd_ds = xr.open_dataset(ctd_netcdf)

ctd_sal = ctd_ds["SAL_ABSOLUTE"].values
ctd_temp = ctd_ds["CONSERVATIVE_TEMP"].values
ctd_castnums = ctd_ds["cast"].values

depth = nitrate_ds["depth"].values
nitrate = nitrate_ds["Nitrate"].values
nitr_lats = nitrate_ds["latitude"].values
nitr_lons = nitrate_ds["longitude"].values
castnums = nitrate_ds["cast"].values


sigma0 = gsw.sigma0(ctd_sal, ctd_temp)


# Extract bathymetry section columns for plotting
section_distances = bathymetry_file['Distance_km']  
bed = bathymetry_file['Bed_Elevation_m'] 
section_lats = bathymetry_file['Latitude']  
section_lons = bathymetry_file['Longitude'] 

start_lon, start_lat = section_lons[0], section_lats[0]
end_lon, end_lat = section_lons.iloc[-1], section_lats.iloc[-1]

# Extract nutrient data
start_row, end_row = 41, 94

nitrate_value = pd.to_numeric(nutrients_file['NO3'].iloc[start_row:end_row])
sample_cast = nutrients_file['St#'].iloc[start_row:end_row]
sample_depth = -pd.to_numeric(nutrients_file['Depth  '].iloc[start_row:end_row])

nitrate_value = nitrate_value.reset_index(drop=True)
sample_cast = sample_cast.reset_index(drop=True)
sample_depth = sample_depth.reset_index(drop=True)



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


# Calculate the distance of each cast along the section and the cast max depth     
    
ctd_distances = []
ctd_maxdepths = []


for i in range(len(castnums)):   

    ctd_distances.append(find_distance(nitr_lats[i], nitr_lons[i], start_lat, start_lon))

    valid = np.where(~np.isnan(ctd_sal[i]))[0]

    if len(valid) > 0:
        ctd_maxdepths.append(float(depth[valid[-1]]))
        
    else: ctd_maxdepths.append(0)

ctd_maxdepths = np.array(ctd_maxdepths)
ctd_distances = np.array(ctd_distances)

start_cast = 11
end_cast = 19
            
# Pad the ctd data so all casts are the same length as nitrate cast data
target_len = nitrate.shape[1]

sal_padded = np.pad(ctd_sal, ((0, 0), (0, target_len - ctd_sal.shape[1])), constant_values=np.nan)
temp_padded = np.pad(ctd_temp, ((0, 0), (0, target_len - ctd_sal.shape[1])), constant_values=np.nan)
sigma0_padded = np.pad(sigma0, ((0, 0), (0, target_len - sigma0.shape[1])), constant_values=np.nan)


# Create the section
idx = np.arange(start_cast-1, end_cast)

section = xr.Dataset(
    data_vars={
        "Absolute Salinity": (["distance", "depth"], sal_padded[idx]),
        "Conservative Temperature": (["distance", "depth"], temp_padded[idx]),
        "Nitrate": (["distance", "depth"], nitrate[idx]),
        "Potential Density Anomaly": (["distance", "depth"], sigma0_padded[idx]),
        "cast_nums": (["distance"], castnums[idx]),
        "max_depth": (["distance"], ctd_maxdepths[idx])
    },
    coords={
        "distance": ctd_distances[idx],
        "depth": depth
    }
)

section = section.sortby("distance")    # Sort by distance


########################################################################################
# Plot salinity section
########################################################################################

nitrate_min = 0
whole_depth_sal = False  ## ENTER WHETHER WHOLE DEPTH OR TOP SEVERAL METERS

if whole_depth_sal == True: 
    text_height = 60
    triangle_height = 40
    nitrate_max = 19.55

else: 
    text_height = 20
    triangle_height = 15
    nitrate_max = 16


# Build meshgrid
X, Y = np.meshgrid(section["distance"], -section["depth"])

levels = np.linspace(nitrate_min, nitrate_max, 200)

fig, ax = plt.subplots(figsize=(12, 6))
cf = ax.contourf(X, Y, section["Nitrate"].T, levels=levels, norm=PowerNorm(gamma=1.3), cmap='YlGn', vmin=nitrate_min, vmax=nitrate_max)
cbar = plt.colorbar(cf, ax=ax, format="%.2f")
cbar.set_label(r"Nitrate [$\mu$M]", labelpad=15)

# Overlay sigma0 contours
sigma_levels = [26, 26.5, 27, 27.5, 28]
smoothed_sigma = savgol_filter(section["Potential Density Anomaly"], window_length=15, polyorder=1)
cs = ax.contour(X, Y, smoothed_sigma.T, levels=sigma_levels, colors='darkslategrey', linewidths=0.7)

label_positions = [(48, -390), (50, -100), (50, -70), (65, -10)]
ax.clabel(cs, inline=True, manual=label_positions, fmt='%.2f', inline_spacing = 20)

plt.tick_params(axis='both', which='major')

for i in range(len(section['distance'])):
    ax.vlines(section['distance'][i], -section['max_depth'][i], triangle_height, zorder = 3, color="k", linewidth = 0.5)
    ax.scatter(section['distance'][i], triangle_height, c="k", marker = "^", s = 50, zorder = 4)
    ax.text(section['distance'][i]-1.5, text_height, str(section['cast_nums'][i].values))
    
    

### Plot nutrient sample values:
    
colormap = plt.colormaps['YlGn']
norm = mcolors.Normalize(vmin=nitrate_min, vmax=nitrate_max)
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)

for j in range(len(sample_cast)-1, -1, -1): #plot shallow samples first

    index = np.where(ctd_castnums==sample_cast[j])[0]
    index = int(index[0])
    # print(sample_cast[j], ctd_castnums[index])
    
    size = nitrate_value[j]*10 + 10
    color = colormap(norm(nitrate_value[j]))
    ax.scatter(ctd_distances[index], sample_depth[j], color=color, edgecolor = 'k', s= size, zorder=4)
    

#################################
    

ax.plot(section_distances, bed, color = "gray")

ax.fill_between(section_distances, bed, min(bed) - 50, color="gray", zorder=2)

if whole_depth_sal == True: 
    ax.text(2, -740, "South")
    ax.text(90, -740, "North")

    ax.set_ylim(-780, 120)
    ax.set_xlim(0, 100)
    

else: 
    ax.set_ylim(-175, 40)
    ax.set_xlim(7, 95)

ax.set_ylabel("Depth [m]")
ax.set_xlabel("Distance Along Section [km]")



