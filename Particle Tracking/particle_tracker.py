#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:32:57 2026

@author: nataliemcgee
"""

import math
from datetime import timedelta
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.cm as cm
import cartopy.feature
import warnings
import parcels
import pyproj

warnings.filterwarnings("ignore", category=UserWarning, module="parcels")

# --- Fieldset ---
filenames = {
    "U": "/Users/nataliemcgee/Desktop/Model Output/Uvel_202106.nc",
    "V": "/Users/nataliemcgee/Desktop/Model Output/Vvel_202106.nc",
    "W": "/Users/nataliemcgee/Desktop/Model Output/Wvel_202106.nc",
}

variables = {
    "U": "Uvel",
    "V": "Vvel",
    "W": "Wvel",
}
dimensions = {"lat": "latitude", "lon": "longitude", 
              "depth": "depths"
              }

start_date = np.datetime64('2021-06-01')

# 1 file containing 30 daily timesteps
timestamps = [np.array([start_date + np.timedelta64(i, 'D') for i in range(30)])]

print(f"Number of files: {len(timestamps)}")           # should be 1
print(f"Timestamps in file: {len(timestamps[0])}")     # should be 30
print(timestamps[0])  

fieldset = parcels.FieldSet.from_netcdf(
    filenames,
    variables,
    dimensions,
    timestamps=timestamps,
    allow_time_extrapolation=True)

fieldset.computeTimeChunk()

fig, ax = plt.subplots(
    figsize=(10, 8),
    subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=-42)}
)

# # Velocity field
# mesh = ax.pcolormesh(
#     fieldset.W.grid.lon,
#     fieldset.W.grid.lat,
#     fieldset.W.data[0, :, :],
#     transform=ccrs.PlateCarree(),
#     zorder=1
# )

# fig.colorbar(mesh, ax=ax, orientation='vertical', pad=0.05, label='U velocity [m/s]')

##### COASTLINE
    
data = np.load("/Users/nataliemcgee/Documents/Upernavik Data/Bathymetry Data/Upernavik_bedmachine_subset1.npz")
x_trim   = data["x"]
y_trim   = data["y"]
bed_trim = data["bed"]

xx_trim, yy_trim = np.meshgrid(x_trim, y_trim)

proj = pyproj.Proj('EPSG:3413')
lon_bath, lat_bath = proj(xx_trim, yy_trim, inverse=True)  # now in lon/lat

cmap = plt.get_cmap('gist_gray').copy()
cmap.set_over('k')

pc = ax.pcolormesh(lon_bath, lat_bath, bed_trim,
                   transform=ccrs.PlateCarree(),
                   cmap=cmap,
                   vmax = -3,
                   zorder = 0)

fig.colorbar(pc, ax=ax, orientation='vertical', pad=0.05, label='Depth')
    
############
depth = 150
pset = parcels.ParticleSet.from_line(
    fieldset=fieldset,
    pclass=parcels.JITParticle,
    size=10,  # releasing 5 particles
    # start=(-54.55, 72.95),  # releasing on a line: the start longitude and latitude
    start=(-54.70, 72.88),
    finish=(-54.45, 73),  # releasing on a line: the end longitude and latitude
    depth = depth,
    time=start_date
)

# Initial particle positions — BEFORE execute()
ax.plot(pset.lon, pset.lat, "ko",
        transform=ccrs.PlateCarree(), zorder=5, label='Start')  

#--- Run simulation ---
output_file = pset.ParticleFile(name="Upernavik_ex.zarr", outputdt=timedelta(hours=1))
pset.execute(
    parcels.AdvectionRK4,
    runtime=timedelta(days=29),
    dt=timedelta(hours=4),
    output_file=output_file,
)

# Final particle positions — AFTER execute()
ax.plot(pset.lon, pset.lat, "rx",
        transform=ccrs.PlateCarree(), zorder=5, label='End')    


# Load and plot full tracks
ds = xr.open_zarr("Upernavik_ex.zarr")
for i in range(ds.dims['trajectory']):
    ax.plot(ds['lon'][i], ds['lat'][i], transform=ccrs.PlateCarree(), linewidth=1, zorder=3)
    

ax.set_extent([-55.15, -54.1, 72.71, 73.05], crs=ccrs.PlateCarree())

ax.text(-55.11, 72.74, f"Depth: {depth} m", transform=ccrs.PlateCarree(), color = "white")

plt.legend(loc="lower left")
plt.show()



