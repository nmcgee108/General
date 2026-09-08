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
from datetime import datetime, timedelta
import pyproj
import xarray as xr

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
    return datetime.strptime(f"{date} {time}", "%Y/%m/%d %H:%M:%S.%f")
    
timestamps = np.array([to_timestamp(d, t) for d, t in zip(dates, times)])
    
start_date, start_time = "2026/07/01", "00:00:00.00"  ## ALL
end_date, end_time   =   "2026/9/05", "00:00:00.00"

# start_date, start_time = "2026/08/03", "22:18:00.00"  ## Northern section
# end_date, end_time   =   "2026/08/05", "00:00:00.00"

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
ax.set_extent([-47.4, -38, 59.4, 60.5], crs=ccrs.PlateCarree())  # full southern end of greenland
#ax.set_extent([-44, -40, 59.5, 60.5], crs=ccrs.PlateCarree()) # Right side of greenland
# ax.set_extent([-48, -44.8, 58.8, 60.5], crs=ccrs.PlateCarree())


# Bathymetry — lon/lat coords, PlateCarree transform
pc = ax.pcolormesh(lon_bath, lat_bath, bed_trim,
                   transform=map_crs,
                   cmap='gist_gray')

# SST scatter — also lon/lat, also PlateCarree

param_to_plot = SSSal

max_val = np.nanmax(param_to_plot)
min_val = np.nanmin(param_to_plot)
#max_val = 11 #temp
min_val = 27 #salinity 
# param_colormap = plt.colormaps['plasma']
param_colormap = plt.colormaps['viridis']
param_norm     = mcolors.Normalize(vmin=min_val, vmax=max_val)
param_sm       = plt.cm.ScalarMappable(cmap=param_colormap, norm=param_norm)
param_sm.set_array([])


valid = np.isfinite(param_to_plot) & np.isfinite(lon) & np.isfinite(lat) & (flow>25)
ax.scatter(lon[valid], lat[valid],
           color=param_colormap(param_norm(param_to_plot[valid])),
           s=40, transform=map_crs, zorder=4)

###### Add ADCP vectors ######

netcdf = "/Users/nataliemcgee/Documents/OSNAP/adcp_data_CF_LS.nc"
ds = xr.open_dataset(netcdf)

depth = 20 #Set the depth of ADCP info we want
print("ADCP Depth:", ds["depth"].sel(depth=depth, method="nearest").values)

adcp_timestamps = ds["timestamp"].values

start_np = np.datetime64(start_timestamp)
end_np   = np.datetime64(end_timestamp)

time_mask = (adcp_timestamps >= start_np) & (adcp_timestamps <= end_np)

adcp_lon = ds["longitude"][time_mask].values
adcp_lat = ds["latitude"][time_mask].values
uvel = ds["Uvel_dt"][time_mask].sel(depth=depth, method="nearest").values
vvel = ds["Vvel_dt"][time_mask].sel(depth=depth, method="nearest").values
adcp_timestamps = adcp_timestamps[time_mask]

# Interpolate data in time, plot every 30 mins
times_interp = np.arange(adcp_timestamps[0], adcp_timestamps[-1], timedelta(minutes = 60))

uvel_interp = np.interp(times_interp.astype('float64'), 
                        adcp_timestamps.astype('float64'), 
                        uvel)

vvel_interp = np.interp(times_interp.astype('float64'), 
                        adcp_timestamps.astype('float64'), 
                        vvel)

lon_interp = np.interp(times_interp.astype('float64'), 
                       adcp_timestamps.astype('float64'), 
                       adcp_lon)

lat_interp = np.interp(times_interp.astype('float64'), 
                       adcp_timestamps.astype('float64'), 
                       adcp_lat)



magnitude = [(u**2+v**2)**(1/2) for u, v in zip(uvel_interp, vvel_interp)]

adcp_colormap = plt.colormaps['Blues']
adcp_norm = mcolors.Normalize(vmin=0.008, vmax=0.7)
adcp_sm = plt.cm.ScalarMappable(cmap=adcp_colormap, norm=adcp_norm)


q=ax.quiver(
    lon_interp,
    lat_interp,
    uvel_interp,
    vvel_interp,
    angles='xy',
    color=adcp_colormap(adcp_norm(magnitude)),
    transform=ccrs.PlateCarree(),
    zorder = 3,
    scale_units='inches', 
    scale=0.3
)

qk1 = ax.quiverkey(
    q, 
    X=0.85, Y=0.75,       # location
    U=0.5,              # vector magnitude represented by the reference arrow
    label='0.5 m/s',
    color=adcp_colormap(0.5),
    labelpos='S',
    labelcolor = "white"      
)

qk2 = ax.quiverkey(
    q, 
    X=0.85, Y=0.9,       # location
    U=0.1,              # vector magnitude represented by the reference arrow
    label='0.1 m/s',
    color=adcp_colormap(0.1),
    labelpos='S'       ,
    labelcolor = "white" 
)

# Colorbars
cbar1 = fig.colorbar(pc, ax=ax, location='bottom', aspect = 30)
cbar1.set_label("Bed Elevation [m]")
cbar2 = fig.colorbar(param_sm, ax=ax, pad=0.02, location='right')
# cbar2.set_label("Temperature [°C]")
cbar2.set_label("Salinity [PSU]")

plt.show()





