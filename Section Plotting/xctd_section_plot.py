#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 11:15:28 2026

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
from matplotlib.lines import Line2D

plt.rcParams['font.size']=14

xctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/XCTD data/adjusted_xctd.nc"
ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"
bathymetry_file = pd.read_csv("/Users/nataliemcgee/Documents/Upernavik Data/Bathymetry Data/northern_trough_bathymetry.csv")

xctd_ds = xr.open_dataset(xctd_netcdf)
ctd_ds = xr.open_dataset(ctd_netcdf)

start_lon, start_lat = -63.1648, 74.5593
end_lon, end_lat = -60.0573, 73.4334

xctd_depth = xctd_ds["depth"].values
xctd_cast = xctd_ds["cast"].values
xctd_sal = xctd_ds["Absolute_Salinity"].values
xctd_temp = xctd_ds["Conservative_Temperature"].values
xctd_lats = xctd_ds["latitude"].values
xctd_lons = xctd_ds["longitude"].values
xctd_castnums = xctd_ds["cast"].values


ctd_depth = ctd_ds["depth"].values
ctd_cast = ctd_ds["cast"].values
ctd_sal = ctd_ds["SAL_ABSOLUTE"].values
ctd_temp = ctd_ds["CONSERVATIVE_TEMP"].values
ctd_lats = ctd_ds["LAT"].values
ctd_lons = ctd_ds["LON"].values
ctd_castnums = ctd_ds["cast"].values


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

xctd_distances = []
xctd_maxdepths = []

for i in range(0, 7):    # iterate through XCTD casts 1-7 
    xctd_distances.append(find_distance(xctd_lats[i], xctd_lons[i], start_lat, start_lon))

    valid = np.where(~np.isnan(xctd_sal[i]))[0]

    if len(valid) > 0:
        xctd_maxdepths.append(float(xctd_depth[valid[-1]]))
        
        
# Calculate the distance of each cast along the section and the cast max depth     
    
ctd_distances = []
ctd_maxdepths = []

for i in range(25, 33):    # iterate through CTD casts 26-33
    ctd_distances.append(find_distance(ctd_lats[i], ctd_lons[i], start_lat, start_lon))

    valid = np.where(~np.isnan(ctd_sal[i]))[0]

    if len(valid) > 0:
        ctd_maxdepths.append(float(ctd_depth[valid[-1]]))

# Merge the xctd and ctd datasets

# Pad the xctd data so all casts are the same length as CTD cast data
target_len = ctd_sal.shape[1]

xctd_sal_padded = np.pad(xctd_sal[0:7], ((0, 0), (0, target_len - xctd_sal.shape[1])), constant_values=np.nan)
xctd_temp_padded = np.pad(xctd_temp[0:7], ((0, 0), (0, target_len - xctd_sal.shape[1])), constant_values=np.nan)

sal_combined = np.concatenate((xctd_sal_padded, ctd_sal[25:33]), axis = 0)
temp_combined = np.concatenate((xctd_temp_padded, ctd_temp[25:33]), axis = 0)
sigma0_combined = gsw.sigma0(sal_combined, temp_combined)


distances_combined = np.concatenate((xctd_distances, ctd_distances))
maxdepths_combined = np.concatenate((xctd_maxdepths, ctd_maxdepths))
castnums = np.concatenate((xctd_castnums[0:7], ctd_castnums[25:33]))
casttypes = ["xctd", "xctd", "xctd", "xctd", "xctd", "xctd", "xctd", 
             "ctd", "ctd", "ctd", "ctd", "ctd", "ctd", "ctd", "ctd", ]

            
section = xr.Dataset(
    data_vars={
        "Absolute Salinity": (["distance", "depth"], sal_combined),
        "Conservative Temperature": (["distance", "depth"], temp_combined),
        "Potential Density Anomaly": (["distance", "depth"], sigma0_combined),
        "cast_nums": (["distance"], castnums),
        "cast_types": (["distance"], casttypes),
        "max_depth": (["distance"], maxdepths_combined)},
    coords={
        "distance": distances_combined,
        "depth": ctd_depth})

section = section.sortby("distance")    # Sort by distance


# Plot bathymetry

# Extract columns for plotting
section_distances = bathymetry_file['Distance_km']  
bed = bathymetry_file['Bed_Elevation_m']


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

label_positions = [(65, -400), (60, -100), (55, -70), (40, -30)]
ax.clabel(cs, inline=True, manual=label_positions, fmt='%.2f', inline_spacing = 20)

plt.tick_params(axis='both', which='major')

for i in range(len(section['distance'])):
    if section["cast_types"][i]=="xctd":
        ax.vlines(section['distance'][i], -section['max_depth'][i], 20, zorder = 3, color="k", linewidth = 0.5, linestyle = ":")
        ax.scatter(section['distance'][i], 20, c="gold", marker = "^", s = 50, zorder = 4)
        
    if section["cast_types"][i]=="ctd":
        ax.vlines(section['distance'][i], -section['max_depth'][i], 20, zorder = 3, color="k", linewidth = 0.5)
        ax.scatter(section['distance'][i], 20, c="red", marker = "o", s = 50, zorder = 4)
    
ax.plot(section_distances, bed, color = "gray")
ax.fill_between(section_distances, bed, min(bed) - 50, color="gray", zorder=0)

custom_handles = [
    Line2D([0], [0], ls='', markeredgecolor='gold', markerfacecolor = "gold", marker = "^", label='XCTD'),
    Line2D([0], [0], ls='', markeredgecolor='red', markerfacecolor = "red",  marker = "o", label='Ship-based CTD')]

ax.legend(handles=custom_handles, loc='lower right')

ax.text(5, -700, "North")

ax.set_ylim(-750, 70)
ax.set_xlim(0, 155)
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
cf = ax.contourf(X, Y, section["Conservative Temperature"].T, levels=levels, norm=PowerNorm(gamma=0.55), cmap='RdBu_r')
cbar = plt.colorbar(cf, ax=ax, format="%.2f")
cbar.set_label("Conservative Temperature [°C]", labelpad=15)

# Overlay sigma0 contours
sigma_levels = [26, 26.5, 27, 27.5, 28]
smoothed_sigma = savgol_filter(section["Potential Density Anomaly"], window_length=15, polyorder=1)
cs = ax.contour(X, Y, smoothed_sigma.T, levels=sigma_levels, colors='white', linewidths=1)

label_positions = [(65, -400), (60, -100), (55, -70), (40, -30)]
ax.clabel(cs, inline=True, manual=label_positions, fmt='%.2f', inline_spacing = 20)

plt.tick_params(axis='both', which='major')

for i in range(len(section['distance'])):
    if section["cast_types"][i]=="xctd":
        ax.vlines(section['distance'][i], -section['max_depth'][i], 20, zorder = 3, color="k", linewidth = 0.5, linestyle = ":")
        ax.scatter(section['distance'][i], 20, c="gold", marker = "^", s = 50, zorder = 4)
        
    if section["cast_types"][i]=="ctd":
        ax.vlines(section['distance'][i], -section['max_depth'][i], 20, zorder = 3, color="k", linewidth = 0.5)
        ax.scatter(section['distance'][i], 20, c="red", marker = "o", s = 50, zorder = 4)
        
ax.text(5, -700, "North")
    
ax.plot(section_distances, bed, color = "gray")

ax.fill_between(section_distances, bed, min(bed) - 50, color="gray", zorder=0)
custom_handles = [
    Line2D([0], [0], ls='', markeredgecolor='gold', markerfacecolor = "gold", marker = "^", label='XCTD'),
    Line2D([0], [0], ls='', markeredgecolor='red', markerfacecolor = "red",  marker = "o", label='Ship-based CTD')]

ax.legend(handles=custom_handles, loc='lower right')

ax.set_ylim(-750, 70)
ax.set_xlim(0, 155)
ax.set_ylabel("Depth [m]")
ax.set_xlabel("Distance Along Section [km]")





    