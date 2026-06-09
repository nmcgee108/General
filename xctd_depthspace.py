#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 12:44:15 2026

@author: nataliemcgee
"""

import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math

ctd_netcdf = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"
xctd_netcdf = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/XCTD data/converted_xctd.nc"

xctd_ds = xr.open_dataset(xctd_netcdf)
ctd_ds = xr.open_dataset(ctd_netcdf)


ctd_depth = ctd_ds["depth"].values
ctd_sal= ctd_ds["SAL_ABSOLUTE"].values
ctd_temp= ctd_ds["CONSERVATIVE_TEMP"].values
ctd_lats = ctd_ds["LAT"].values
ctd_lons = -ctd_ds["LON"].values
ctd_nums = ctd_ds["CAST_NUM"].values


xctd_depth = xctd_ds["depth"].values
xctd_sal = xctd_ds["Salinity"].values
xctd_temp = xctd_ds["Temperature"].values
xctd_lats = xctd_ds["latitude"].values
xctd_lons = xctd_ds["longitude"].values

fig, axes = plt.subplots(1,2,figsize=(9, 8), sharey=True)

def find_distance(lat1, lon1, lat2, lon2):
    
    # lat and lon of station 4:
    lat1 = lat1
    lon1 = lon1
    
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


xctd_cast = 7


colormap = plt.colormaps['viridis']
norm = mcolors.Normalize(vmin=0, vmax=60)
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)


print("Closest Casts and Distance from XCTD [km]:")

for i in range(25,56):
    
    dist = find_distance(xctd_lats[xctd_cast-1], xctd_lons[xctd_cast-1], ctd_lats[i], ctd_lons[i])
    color = colormap(norm(dist))
    
    if dist > 60:     # only plot the casts within 50 km
        continue

    print(ctd_nums[i],",",dist)

    axes[0].plot(ctd_sal[i], -ctd_depth, color = color, label = ctd_nums[i])
    axes[0].set_xlim(32, 35)
    
    axes[1].plot(ctd_temp[i], -ctd_depth, color = color)
    axes[1].set_xlim(-2, 4)
    
    
    
axes[0].plot(xctd_sal.T[xctd_cast-1]+0.2, -xctd_depth, color = "red")
axes[1].plot(xctd_temp.T[xctd_cast-1], -xctd_depth, color = "red")

cbar_ax = fig.add_axes([0.95,0.1,0.05,0.77])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Distance from XCTD [km]', fontsize = 14)

axes[0].set_ylabel("Depth [m]")
axes[0].set_xlabel("SA [g/kg]")
axes[1].set_xlabel("CT [°C]")

axes[0].legend(title = "Closest Casts")

fig.suptitle(f"XCTD Cast {xctd_cast} and Closest CTD casts")























