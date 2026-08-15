#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 02:35:06 2026

@author: nataliemcgee
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import glob
import pandas as pd
import matplotlib.colors as mcolors
from datetime import datetime
import pyproj

plt.rcParams['font.size'] = 18

# ------------------------------------------------------------------ #
#  RULE: everything is in lon/lat, transform=PlateCarree() everywhere #
# ------------------------------------------------------------------ #

map_crs = ccrs.PlateCarree()

# --- Underway ship data ---
file_pattern = "/Volumes/data_on_memory/underway/proc/AR260*.csv"
file_list = sorted(glob.glob(file_pattern))
print(f"Found {len(file_list)} files")

df = pd.concat([pd.read_csv(f, skiprows=1) for f in file_list], ignore_index=True)
dates = df["DATE_GMT"].values
times = df[" TIME_GMT"].values
lat = df[" Dec_LAT"].values
lon = df[" Dec_LON"].values
multibeam_depth = df[" EM124"].values #? [m]
SSSal = df[" SBE45S"].values# PSU
SSTemp = df[" SBE48T"].values
flow = df[" FLOW"]


def to_timestamp(date, time):
    timestamp = datetime.strptime(f"{date} {time}", "%Y/%m/%d %H:%M:%S.%f").timestamp()
    return timestamp
    
timestamps = np.array([to_timestamp(d, t) for d, t in zip(dates, times)])
    
    
start_date, start_time = "2026/08/03", "22:18:00.00"  ## Northern section
end_date, end_time   =   "2026/08/05", "00:00:00.00"

# start_date, start_time = "2026/08/11", "00:30:00.00" ## PCS Section 1
# end_date, end_time   =   "2026/08/11", "18:00:00.00"

# start_date, start_time = "2026/08/07", "21:30:00.00" ## PCS Section 2
# end_date, end_time   =   "2026/08/10", "23:57:00.00"


start_timestamp = to_timestamp(start_date, start_time)
end_timestamp = to_timestamp(end_date, end_time)

time_mask = (timestamps >= start_timestamp) & (timestamps <= end_timestamp)

lat = lat[time_mask]
lon = lon[time_mask]
multibeam_depth = multibeam_depth[time_mask]
SSSal = SSSal[time_mask]
SSTemp = SSTemp[time_mask]
flow = flow[time_mask]


# --- Bathymetry: convert EPSG:3413 x/y → lon/lat ONCE, then forget about it ---
data = np.load("/Users/nataliemcgee/Documents/Upernavik Data/Bathymetry Data/OSNAP_bedmachine_subset1.npz")
x_trim   = data["x"]
y_trim   = data["y"]
bed_trim = data["bed"]

xx_trim, yy_trim = np.meshgrid(x_trim, y_trim)

proj = pyproj.Proj('EPSG:3413')
lon_bath, lat_bath = proj(xx_trim, yy_trim, inverse=True)  # now in lon/lat


# --- Plot ---
fig = plt.figure(figsize=(23, 5), layout='compressed')

# Change ONLY this line to change how the map looks
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-42))

ax.coastlines(resolution='10m')
#ax.set_extent([-46.2, -40, 59.52, 60.5], crs=ccrs.PlateCarree())  # extent always in lon/lat
ax.set_extent([-44, -40, 59.5, 60.5], crs=ccrs.PlateCarree())


# Bathymetry — lon/lat coords, PlateCarree transform
pc = ax.pcolormesh(lon_bath, lat_bath, bed_trim,
                   transform=map_crs,
                   cmap='gist_gray')

# SST scatter — also lon/lat, also PlateCarree

param_to_plot = SSTemp

max_val = np.nanmax(param_to_plot)
min_val = np.nanmin(param_to_plot)
max_val = 11
param_colormap = plt.colormaps['plasma']
param_norm     = mcolors.Normalize(vmin=min_val, vmax=max_val)
param_sm       = plt.cm.ScalarMappable(cmap=param_colormap, norm=param_norm)
param_sm.set_array([])


valid = np.isfinite(param_to_plot) & np.isfinite(lon) & np.isfinite(lat) & (flow>25)
ax.scatter(lon[valid], lat[valid],
           color=param_colormap(param_norm(param_to_plot[valid])),
           s=40, transform=map_crs, zorder=4)

# Colorbars
cbar1 = fig.colorbar(pc, ax=ax, location='bottom', aspect = 30)
cbar1.set_label("Bed Elevation [m]")
cbar2 = fig.colorbar(param_sm, ax=ax, pad=0.02, location='right')
cbar2.set_label("Temperature [°C]")

plt.show()





