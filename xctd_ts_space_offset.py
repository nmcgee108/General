#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 11:15:23 2026

@author: nataliemcgee
"""


import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math
import numpy as np
import gsw

ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"
xctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/XCTD data/converted_xctd.nc"

xctd_ds = xr.open_dataset(xctd_netcdf)
ctd_ds = xr.open_dataset(ctd_netcdf)


ctd_depth = ctd_ds["depth"].values
ctd_sal= ctd_ds["SAL_ABSOLUTE"].values
ctd_temp= ctd_ds["CONSERVATIVE_TEMP"].values
ctd_lats = ctd_ds["LAT"].values
ctd_lons = -ctd_ds["LON"].values
ctd_nums = ctd_ds["CAST_NUM"].values


xctd_depth = xctd_ds["depth"].values
xctd_pract_sal = xctd_ds["Salinity"].values
xctd_temp = xctd_ds["Temperature"].values
xctd_lats = xctd_ds["latitude"].values
xctd_lons = xctd_ds["longitude"].values

depth2d = np.tile(xctd_depth[:, np.newaxis], (1, len(xctd_lats)))
xctd_pres = gsw.p_from_z(-depth2d, xctd_lats)
xctd_sal = gsw.SA_from_SP(xctd_pract_sal, xctd_pres, xctd_lons, xctd_lats)
print('done')


fig, axes = plt.subplots(1,1,figsize=(7, 6), sharey=True)

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
dist_cutoff = 40
offset = 0.05


colormap = plt.colormaps['viridis']
norm = mcolors.Normalize(vmin=0, vmax=dist_cutoff)
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)

sx = np.arange(30, 36, 0.5)
ty = np.arange(-2, 5, 0.5)

# Create grid of salinity and temperature
S, T = np.meshgrid(sx, ty)
PDEN = gsw.rho(S, T, 0) - 1000  # Potential density anomaly at reference pressure = 0

# Plot the isopycnal contour lines
contour = axes.contour(S, T, PDEN, levels=[25.5, 26, 26.5, 27, 27.5, 28], colors='grey')

# Manually label contours
label_positions = [(32.8, 2.1), (33.1, 2.1), (33.8, 2.1), (34.5, 1.4), (35.0, 1)]
plt.clabel(contour, inline=True, manual=label_positions, fmt='%1.1f')

print("Closest Casts and Distance from XCTD [km]:")

close_casts = []

for i in range(25,56):
    
    dist = find_distance(xctd_lats[xctd_cast-1], xctd_lons[xctd_cast-1], ctd_lats[i], ctd_lons[i])
    color = colormap(norm(dist))
    
    if dist > dist_cutoff:     # only plot the casts within 50 km
        continue

    print(ctd_nums[i],",",dist)
    close_casts.append(i)

    axes.scatter(ctd_sal[i], ctd_temp[i], color = color, label = ctd_nums[i], s = 4)
    axes.set_xlim(32, 35)



# # STEP 1: Identify the XCTD and CTD Deep Water (adjust bounds as necessary)
# xctd_sal_deep, xctd_temp_deep = xctd_sal.T[xctd_cast-1][-500:],  xctd_temp.T[xctd_cast-1][-500:]
 
axes.scatter(xctd_sal.T[xctd_cast-1]+offset, xctd_temp.T[xctd_cast-1], color = "red", s = 4)
# axes.scatter(xctd_sal_deep, xctd_temp_deep, color = "red", s = 4)

# # STEP 2: Calculate a reference temperature by taking the deep water average
# xctd_dw_temp = np.nanmean(xctd_temp_deep)
# axes.hlines(xctd_dw_temp, 33.5, 34.75, linestyle="dashed", color = "k")

# print( )
# print("Deep Water Reference Temp:", xctd_dw_temp)

# # STEP 3: Calculate the average salinity of the CTD water near this temperature
# for i in close_casts:
#     pass

    
    
cbar_ax = fig.add_axes([0.95,0.1,0.05,0.77])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Distance from XCTD [km]')

axes.set_ylabel("CT [°C]")
axes.set_xlabel("SA [g/kg]")

axes.set_ylim(1, 3)
axes.set_xlim(34,34.8)


axes.legend(title = "Closest Casts")

fig.suptitle(f"XCTD Cast {xctd_cast} and Closest CTD casts"+"\nSA Offset = +"+str(offset))










