#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 05:55:46 2026

@author: nataliemcgee
"""

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import glob
import pandas as pd
import matplotlib.colors as mcolors
import pyproj

plt.rcParams['font.size'] = 18

# --- Underway data ---
file_pattern = "/Volumes/data_on_memory/underway/proc/AR260*.csv"
file_list = sorted(glob.glob(file_pattern))
df = pd.concat([pd.read_csv(f, skiprows=1) for f in file_list], ignore_index=True)

lat    = df[" Dec_LAT"].values
lon    = df[" Dec_LON"].values
SSTemp = df[" SBE48T"].values

# --- Bathymetry ---
ds  = xr.open_dataset("/Users/nataliemcgee/Desktop/Large Datasets/BedMachineGreenland-v5-002.nc")
bed = ds['bed'].values
x   = ds['x'].values    # 1D
y   = ds['y'].values    # 1D

proj = pyproj.Proj('EPSG:3413')

# --- KEY SIMPLIFICATION ---
# Instead of converting the whole grid to lon/lat,
# convert your desired lon/lat box corners into native x/y and slice directly

lon_min, lon_max = -63.5, -55
lat_min, lat_max =  72.6,  75.5

# Convert all 4 corners to native x/y to get a safe bounding box
cx, cy = proj([lon_min, lon_max, lon_min, lon_max],
              [lat_min, lat_min, lat_max, lat_max])

x_trim   = x[(x >= min(cx)) & (x <= max(cx))]          # slice 1D arrays directly
y_trim   = y[(y >= min(cy)) & (y <= max(cy))]
x_mask   = (x >= min(cx)) & (x <= max(cx))
y_mask   = (y >= min(cy)) & (y <= max(cy))
bed_trim = bed[np.ix_(y_mask, x_mask)]                  # slice 2D bed array

# Convert trimmed grid to lon/lat for plotting
xx_trim, yy_trim = np.meshgrid(x_trim, y_trim)
lon_bath, lat_bath = proj(xx_trim, yy_trim, inverse=True)

# --- SAVE SUBSET ---
np.savez(
    "Upernavik_bedmachine_subset1.npz",
    bed=bed_trim,
    x=x_trim,
    y=y_trim
)
print(f"Saved subset: {bed_trim.shape}")

# --- Plot ---
fig = plt.figure(figsize=(19, 6))
ax  = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-42))

pc = ax.pcolormesh(lon_bath, lat_bath, bed_trim,
                   transform=ccrs.PlateCarree(),
                   cmap='gist_gray')

ax.coastlines(resolution='10m')
ax.set_extent([-50, -25, 58, 64.5], crs=ccrs.PlateCarree())

valid = np.isfinite(lon) & np.isfinite(lat)
ax.scatter(lon[valid], lat[valid],
           color='red', s=10,
           transform=ccrs.PlateCarree(), zorder=4)

fig.colorbar(pc, ax=ax, pad=0.02, location='right', label='Bed Elevation [m]')
plt.show()



plt.show()