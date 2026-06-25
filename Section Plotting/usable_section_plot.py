#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 14:40:14 2026

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

plt.rcParams['font.size']=14

ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"
bathymetry_file = pd.read_csv("/Users/nataliemcgee/Documents/Upernavik Data/Bathymetry Data/new_fjord_bathymetry1.csv")

ctd_ds = xr.open_dataset(ctd_netcdf)


# Extract bathymetry section columns for plotting
section_distances = bathymetry_file['Distance_km']  
bed = bathymetry_file['Bed_Elevation_m'] 
section_lats = bathymetry_file['Latitude']  
section_lons = bathymetry_file['Longitude'] 

start_lon, start_lat = section_lons[0], section_lats[0]
end_lon, end_lat = section_lons.iloc[-1], section_lats.iloc[-1]

ctd_depth = ctd_ds["depth"].values
ctd_sal = ctd_ds["SAL_ABSOLUTE"].values
ctd_temp = ctd_ds["CONSERVATIVE_TEMP"].values
ctd_lats = ctd_ds["LAT"].values
ctd_lons = ctd_ds["LON"].values
ctd_castnums = ctd_ds["cast"].values +1

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


# Calculate the distance of each cast along the section and the cast max depth     
    
ctd_distances = []
ctd_maxdepths = []


for i in range(len(ctd_castnums)):   

    ctd_distances.append(find_distance(ctd_lats[i], ctd_lons[i], start_lat, start_lon))

    valid = np.where(~np.isnan(ctd_sal[i]))[0]

    if len(valid) > 0:
        ctd_maxdepths.append(float(ctd_depth[valid[-1]]))
        
    else: ctd_maxdepths.append(0)


start_cast = 2
end_cast = 9
            
section = xr.Dataset(
    data_vars={
        "Absolute Salinity": (["distance", "depth"], ctd_sal[start_cast-1:end_cast]),
        "Conservative Temperature": (["distance", "depth"], ctd_temp[start_cast-1:end_cast]),
        "Potential Density Anomaly": (["distance", "depth"], sigma0[start_cast-1:end_cast]),
        "cast_nums": (["distance"], ctd_castnums[start_cast-1:end_cast]),
        "max_depth": (["distance"], ctd_maxdepths[start_cast-1:end_cast])},
    coords={
        "distance": ctd_distances[start_cast-1:end_cast],
        "depth": ctd_depth})

section = section.sortby("distance")    # Sort by distance


######################
# Plot salinity section

SA_min = 32
SA_max = 35

# Build meshgrid
X, Y = np.meshgrid(section["distance"], -section["depth"])

levels = np.linspace(32, 35, 200)

fig, ax = plt.subplots(figsize=(12, 6))
cf = ax.contourf(X, Y, section["Absolute Salinity"].T, levels=levels, norm=PowerNorm(gamma=1), cmap='PuBuGn', vmin=SA_min, vmax=SA_max)
cbar = plt.colorbar(cf, ax=ax, format="%.2f")
cbar.set_label("Absolute Salinity [g/kg]", labelpad=15)

# Overlay sigma0 contours
sigma_levels = [26, 26.5, 27, 27.5, 28]
smoothed_sigma = savgol_filter(section["Potential Density Anomaly"], window_length=15, polyorder=1)
cs = ax.contour(X, Y, smoothed_sigma.T, levels=sigma_levels, colors='white', linewidths=1)

label_positions = [(70, -400), (70, -100), (40, -70), (18, -30)]
ax.clabel(cs, inline=True, manual=label_positions, fmt='%.2f', inline_spacing = 20)

plt.tick_params(axis='both', which='major')

for i in range(len(section['distance'])):

    ax.vlines(section['distance'][i], -section['max_depth'][i], 20, zorder = 3, color="k", linewidth = 0.5)
    ax.scatter(section['distance'][i], 20, c="red", marker = "o", s = 50, zorder = 4)
    ax.text(section['distance'][i]+1, 40, str(section['cast_nums'][i].values))


ax.plot(section_distances, bed, color = "gray")
ax.fill_between(section_distances, bed, min(bed) - 50, color="gray", zorder=0)

ax.text(12, -970, "Glacier")
ax.text(98, -970, "Shelf")

ax.set_ylim(-1000, 100)
ax.set_xlim(0, 100)
ax.invert_xaxis()
ax.set_ylabel("Depth [m]")
ax.set_xlabel("Distance Along Section [km]")

######################
# Plot temp section

T_min = -2
T_max = 4

# Build meshgrid
X, Y = np.meshgrid(section["distance"], -section["depth"])

levels = np.linspace(-2, 5, 200)

fig, ax = plt.subplots(figsize=(12, 6))
cf = ax.contourf(X, Y, section["Conservative Temperature"].T, levels=levels, norm=PowerNorm(gamma=1), cmap='RdBu_r')
cbar = plt.colorbar(cf, ax=ax, format="%.2f")
cbar.set_label("Conservative Temperature [°C]", labelpad=15)

# Overlay sigma0 contours
sigma_levels = [26, 26.5, 27, 27.5, 28]
smoothed_sigma = savgol_filter(section["Potential Density Anomaly"], window_length=15, polyorder=1)
cs = ax.contour(X, Y, smoothed_sigma.T, levels=sigma_levels, colors='dimgray', linewidths=1)

label_positions = [(70, -400), (70, -100), (40, -70), (18, -30)]
ax.clabel(cs, inline=True, manual=label_positions, fmt='%.2f', inline_spacing = 20)

plt.tick_params(axis='both', which='major')

for i in range(len(section['distance'])):

    ax.vlines(section['distance'][i], -section['max_depth'][i], 20, zorder = 3, color="k", linewidth = 0.5)
    ax.scatter(section['distance'][i], 20, c="red", marker = "o", s = 50, zorder = 4)
    ax.text(section['distance'][i]+1, 40, str(section['cast_nums'][i].values))
    
    
ax.plot(section_distances, bed, color = "gray")

ax.fill_between(section_distances, bed, min(bed) - 50, color="gray", zorder=0)

ax.text(12, -970, "Glacier")
ax.text(98, -970, "Shelf")

ax.set_ylim(-1000, 100)
ax.set_xlim(0, 100)
ax.invert_xaxis()
ax.set_ylabel("Depth [m]")
ax.set_xlabel("Distance Along Section [km]")





    