#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 11:42:32 2026

@author: nataliemcgee
"""

import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import gsw
import numpy as np
from netCDF4 import Dataset

plt.figure(figsize=(7, 7))
plt.rcParams['font.size'] = 13

argo_netcdf = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/argo_float_data/GL_PR_PF_2904018.nc"

argo_ds = xr.open_dataset(argo_netcdf)

#print(argo_ds["POSITION_QC"])
#print(list(argo_ds.variables))


month_bounds = [["2024-03-01", "2024-03-31"],
                ["2024-04-01", "2024-04-30"],
                ["2024-05-01", "2024-05-31"], 
                ["2024-06-01", "2024-06-30"],
                ["2024-07-01", "2024-07-31"],
                ["2024-08-01", "2024-08-31"]]


colormap = plt.colormaps['viridis'] # Choose colormap
norm = mcolors.Normalize(vmin=0, vmax=5)  # Normalize the colormap with a max and min value
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm) 

labels = ["March", "April", "May", "June", "July", "August"]

argo_lats = []
argo_lons = []

which_cast = 2


for i in range(6):
    month_ds = argo_ds.sel(TIME=slice(month_bounds[i][0], month_bounds[i][1]))

    temp = month_ds["TEMP_ADJUSTED"].isel(TIME = which_cast)
    pres = month_ds["PRES_ADJUSTED"].isel(TIME = which_cast)
    psal = month_ds["PSAL_ADJUSTED"].isel(TIME = which_cast)
    
    temp_qc = month_ds["TEMP_ADJUSTED_QC"].isel(TIME = which_cast)
    pres_qc = month_ds["PRES_ADJUSTED_QC"].isel(TIME = which_cast)
    psal_qc = month_ds["PSAL_ADJUSTED_QC"].isel(TIME = which_cast)
    
    good = np.isin(temp_qc, [1, 2]) & \
           np.isin(psal_qc, [1, 2]) & \
           np.isin(pres_qc, [1, 2]) # only plot points that are "good" or "probably good"
    
    
    lon = month_ds["LONGITUDE"].isel(TIME=which_cast).item()
    lat = month_ds["LATITUDE"].isel(TIME=which_cast).item()
    SA = gsw.SA_from_SP(psal, pres, lon, lat)
    
    argo_lats.append(lat)
    argo_lons.append(lon)

    color = colormap(norm(i))

    plt.scatter(SA.where(good), temp.where(good), color=color, s=1.5, label = labels[i])
    
# Plot the isopycnal contour lines
minS = 32
maxS = 35
minT = -2
maxT = 4

sx = np.arange(minS, maxS + 0.1, 0.1)
ty = np.arange(minT, maxT + 0.1, 0.1)

# Create grid of salinity and temperature
S, T = np.meshgrid(sx, ty)
PDEN = gsw.rho(S, T, 0) - 1000  # Potential density anomaly at reference pressure = 0

contour = plt.contour(S, T, PDEN, levels=[25.5, 26, 26.5, 27, 27.5, 28], colors='grey')

# Manually label contours
label_positions = [(32.8, 3.1), (33.1, 2.5), (33.8, 2.9), (34.5, 1.4), (35.0, 1)]
plt.clabel(contour, inline=True, manual=label_positions, fmt='%1.1f')

    
# Add cast 9 from 2024 shipboard CTD casts for reference

uc_patch_file = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Final CTD Datasets/uc_patch_dataset_new.nc"

ncfile = Dataset(uc_patch_file, 'r')

for i in [2, 8]:
    #NOTE: indexes are as such: CT[7] is cast 8
    SA_shipboard = ncfile.variables['Absolute Salinity'][i]
    CT_shipboard = ncfile.variables['Conservative Temp'][i]
    
    if i == 2:
        shipboard_color = "orange"
        shipboard_label = "Fjord Head Cast 3 (Jul. 31)"
    
    if i == 8:
        shipboard_color = "crimson"
        shipboard_label = "Fjord Mouth Cast 9 (Aug. 1)"
        
    plt.plot(SA_shipboard, CT_shipboard, color = shipboard_color, label = shipboard_label)
    
ncfile.close()

plt.xlabel("Absolute Salinity [g/kg]")
plt.ylabel("Adjusted Temperature [°C]")

plt.xlim(33,35)
plt.ylim(-2,4)
plt.legend()

################################################################################
#  Plot map next
################################################################################


import cartopy.crs as ccrs
from scipy.signal import savgol_filter

data = np.load("/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Bathymetry Data/upernavik_troughandfjord_bedmachine_subset.npz")

bed_trim = data["bed"]
lat_trim = data["lat"]
lon_trim = data["lon"]

fig = plt.figure(figsize=(15, 6))
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-55))

# Plot using PlateCarree (lon/lat) coords
bed_trim_smooth = savgol_filter(bed_trim, window_length=7, polyorder=2, axis=0)
pc = ax.pcolormesh(lon_trim, lat_trim, bed_trim, transform=ccrs.PlateCarree(), cmap='gist_gray')


ax.coastlines() # draw coastlines


ax.set_extent([-61, -54, 73.8, 72.65], crs=ccrs.PlateCarree()) # Set map limits in lon/lat


# Add colorbar 
cbar = fig.colorbar(pc, ax=ax, pad=0.02)
cbar.set_label("Bed Elevation [m]")

#fjord
ctd_latitudes = [
    72.91,
    72.8849,
    72.88805,
    72.8822,
    72.91523,
    72.9404,
    73.00466,
    73.017833,
    73.127,
    72.7735,   # first shelf/trough cast
    72.83995,
    72.91643,
    72.9867,
    73.05838,
    73.12492,
    73.19801,
    73.27043,
    73.3421]

ctd_longitudes = [
    -55.01,
    -55.0707,
    -55.0057,
    -54.9242,
    -55.454,
    -55.6726,
    -56.08366,
    -56.330163,
    -57.3294,
    -58.58975,  # first shelf/trough cast 
    -58.7702,
    -58.92601,
    -59.146,
    -59.34056,
    -59.504433,
    -59.70404,
    -59.90671,
    -60.076983]


ax.scatter(ctd_longitudes, ctd_latitudes, color='k', s=25, transform=ccrs.PlateCarree(), zorder = 4, label="Shipboard CTD Stations (August)")
ax.scatter(-57.3294, 73.127, color='crimson', s=170, marker = "*", transform=ccrs.PlateCarree(), zorder = 3)
ax.scatter(-54.9242, 72.8822, color='orange', s=170, marker = "*", transform=ccrs.PlateCarree(), zorder = 3)

#annotate shipboard casts
ax.annotate("Cast 9 (Aug. 1)", 
            xy=(-57.3294, 73.127), 
            xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
            xytext=(-57.4, 72.8),
            textcoords=ccrs.PlateCarree()._as_mpl_transform(ax),
            color = "crimson",
            arrowprops=dict(arrowstyle="->"),
            fontsize=13,
            zorder=6)

ax.annotate("Cast 3 (Jul. 31)", 
            xy=(-54.9242, 72.8822), 
            xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
            xytext=(-55.25, 73.1),
            textcoords=ccrs.PlateCarree()._as_mpl_transform(ax),
            color = "orange",
            arrowprops=dict(arrowstyle="->"),
            fontsize=13,
            zorder=6)


offsets = [
    (-0.7, 0.1),
    (-0.3, 0.2),
    (0.15, 0.2),
    (0.35, 0.1),
    (0.5, 0),
    (0.2, -0.25)]

for i in range(6):
    color = colormap(norm(i))

    ax.scatter(
        argo_lons[i],
        argo_lats[i],
        color=color,
        s=40,
        transform=ccrs.PlateCarree(),
        zorder=5)
    
    dx, dy = offsets[i]

    ax.annotate(
        labels[i],
        xy=(argo_lons[i], argo_lats[i]),
        xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),

        xytext=(argo_lons[i] + dx, argo_lats[i] + dy),
        textcoords=ccrs.PlateCarree()._as_mpl_transform(ax),
        color = color,
        arrowprops=dict(arrowstyle="->"),
        fontsize=13,
        zorder=6)
    
    


plt.legend(loc = 'upper left')
plt.show()

