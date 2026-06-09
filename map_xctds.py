#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 12:14:08 2026

@author: nataliemcgee
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import cartopy.crs as ccrs


ctd_netcdf = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"
xctd_netcdf = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/XCTD data/converted_xctd.nc"

xctd_ds = xr.open_dataset(xctd_netcdf)
ctd_ds = xr.open_dataset(ctd_netcdf)

# print(xctd_ds["latitude"].values)
# print(xctd_ds["longitude"].values)


data = np.load("/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Bathymetry Data/melvillebay_bedmachine_subset1.npz")

bed_trim = data["bed"]
lat_trim = data["lat"]
lon_trim = data["lon"]

fig = plt.figure(figsize=(15, 10))
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-55))

# Plot using PlateCarree (lon/lat) coords
bed_trim_smooth = savgol_filter(bed_trim, window_length=7, polyorder=2, axis=0)
pc = ax.pcolormesh(lon_trim, lat_trim, bed_trim, transform=ccrs.PlateCarree(), cmap='gist_gray')


ax.coastlines() # draw coastlines


ax.set_extent([-63, -56.5, 75.2, 73.2], crs=ccrs.PlateCarree()) # Set map limits in lon/lat


# Add colorbar 
cbar = fig.colorbar(pc, ax=ax, pad=0.02)
cbar.set_label("Bed Elevation [m]")

#xctds

x_lats = xctd_ds["latitude"].values
x_lons = xctd_ds["longitude"].values

#fjord
ctd_latitudes = ctd_ds["LAT"].values
ctd_longitudes = -ctd_ds["LON"].values

print(ctd_latitudes)
print(ctd_longitudes)

ax.scatter(ctd_longitudes, ctd_latitudes, color='k', s=25, transform=ccrs.PlateCarree(), zorder = 4, label="Shipboard CTD Stations (August)")
ax.scatter(x_lons, x_lats, color='red', s=25, transform=ccrs.PlateCarree(), zorder = 4, label="XCTD Stations")


plt.legend(loc = 'upper left')
plt.show()

