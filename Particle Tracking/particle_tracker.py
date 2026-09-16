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
import time
t0 = time.time()

warnings.filterwarnings("ignore", category=UserWarning, module="parcels")

# --- Fieldset ---
filenames = {
    "U": "/Users/nataliemcgee/Desktop/Model Output/Uvel_202108.nc",
    "V": "/Users/nataliemcgee/Desktop/Model Output/Vvel_202108.nc",
    "W": "/Users/nataliemcgee/Desktop/Model Output/Wvel_202108.nc",
}

variables = {
    "U": "Uvel",
    "V": "Vvel",
    "W": "Wvel",
}
dimensions = {"lat": "latitude", "lon": "longitude", 
              "depth": "depths"
              }

#start_date = np.datetime64('2021-06-01')
# start_date = np.datetime64('2021-07-01')
start_date = np.datetime64('2021-08-01')

# 1 file containing 30 daily timesteps
timestamps = [np.array([start_date + np.timedelta64(i, 'D') for i in range(31)])]

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

print("fieldset ready:", time.time() - t0)

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

#PARTICLE DEPTH
depth = 3

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
                   vmax = -depth,
                   zorder = 0)

print("bathymetry plotted:", time.time() - t0)

fig.colorbar(pc, ax=ax, orientation='vertical', pad=0.05, label='Depth')
    
############
n_horiz = 5    # particles along the line
n_depth = 3     # particles in depth

# lons = np.linspace(-54.53, -54.427, n_horiz)
# lats = np.linspace(72.98, 72.967, n_horiz) # GLACIER N

# lons = np.linspace(-54.43, -54.40, n_horiz)
# lats = np.linspace(72.926, 72.955, n_horiz) # GLACIER C

lons = np.linspace(-54.55, -54.49, n_horiz)
lats = np.linspace(72.850, 72.89, n_horiz) # GLACIER S

# the depths
depths = np.linspace(1, 200, n_depth)
print(depths)

# build the full grid: every horizontal point repeated at every depth
lon_grid = np.tile(lons, n_depth)      # [l0,l1,l2,l3,l4, l0,l1,l2,l3,l4, ...]
lat_grid = np.tile(lats, n_depth)
depth_grid = np.repeat(depths, n_horiz)  # [100,100,100,100,100, 350,350,..., 600,...]

pset = parcels.ParticleSet(
    fieldset=fieldset,
    pclass=parcels.JITParticle,
    lon=lon_grid,
    lat=lat_grid,
    depth=depth_grid,
    time=start_date,   # scalar time broadcasts to all particles
    repeatdt=timedelta(days=5)
)


# Initial particle positions — BEFORE execute()
ax.plot(pset.lon, pset.lat, "ko",
        transform=ccrs.PlateCarree(), zorder=5, label='Start')  

# # --- Run simulation ---
# output_file = pset.ParticleFile(name="/Users/nataliemcgee/Desktop/Model Output/Upernavik_ex.zarr", outputdt=timedelta(hours=1))
# pset.execute(
#     parcels.AdvectionRK4,
#     runtime=timedelta(days=29),
#     dt=timedelta(hours=1),
#     output_file=output_file,
# )

# Final particle positions — AFTER execute()
# ax.plot(pset.lon, pset.lat, "rx",
#         transform=ccrs.PlateCarree(), zorder=5, label='End')    


# Load and plot full tracks
ds = xr.open_zarr("/Users/nataliemcgee/Desktop/Model Output/Upernavik_ex.zarr")

print("zarr opened:", time.time() - t0)


lon = ds['lon'].values   # shape (n_traj, n_obs), loaded into memory once
lat = ds['lat'].values

print("lon/lat loaded, n_traj =", lon.shape[0], "time:", time.time() - t0)

def find_distance_vectorized(lon1, lat1, lon2, lat2):
    R = 6371.0  # km
    phi1 = np.radians(lat1)
    phi2 = np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)

    a = np.sin(dphi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(dlambda/2)**2
    c = 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))
    return R * c


# glacier outlet location
lon1_fixed, lat1_fixed = -54.3532829284668, 72.85934448242188 # Glacier S
# lon1_fixed, lat1_fixed = -54.447349548339844, 72.9875717163086 # Glacier N


dist = find_distance_vectorized(lon, lat, lon1_fixed, lat1_fixed)  # shape (n_traj, n_obs)
tags = np.any(dist < 1, axis=1)   # True if particle ever gets within 1 km

print(tags)

for i in range(lon.shape[0]):
    if tags[i]==True:
        ax.plot(lon[i], lat[i], transform=ccrs.PlateCarree(), color = "green",
                linewidth=1, zorder=3)
    else: ax.plot(lon[i], lat[i], transform=ccrs.PlateCarree(), color = "orange",
                  linewidth=1, zorder=3)
    
print("tracks plotted:", time.time() - t0)


ax.plot(-54.447349548339844, 72.9875717163086, "bo",
        transform=ccrs.PlateCarree(), zorder=5)    

ax.plot(-54.36542892456055, 72.93659210205078, "bo",
        transform=ccrs.PlateCarree(), zorder=5)    

ax.plot(-54.3532829284668, 72.85934448242188, "bo",
        transform=ccrs.PlateCarree(), zorder=5, label='Outlets')  

# ax.plot([-54.43, -54.40], [72.926, 72.955], transform=ccrs.PlateCarree(), zorder = 5, color = "red")#glacier C
# ax.plot([-54.53, -54.427], [72.98, 72.967], transform=ccrs.PlateCarree(), zorder = 5, color = "red") #glaicer N
ax.plot([-54.55, -54.49], [72.850, 72.89], transform=ccrs.PlateCarree(), zorder = 5, color = "red") #glaicer N


ax.set_extent([-55.15, -54.1, 72.71, 73.05], crs=ccrs.PlateCarree())

ax.text(-55.11, 72.74, f"Depth: {depth} m", transform=ccrs.PlateCarree(), color = "white")

plt.legend(loc="lower left")
plt.show()



