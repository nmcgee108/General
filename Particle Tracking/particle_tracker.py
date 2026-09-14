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

fig.colorbar(pc, ax=ax, orientation='vertical', pad=0.05, label='Depth')
    
############
n_horiz = 5    # particles along the line
n_depth = 5     # particles in depth

# the line in lon/lat
# lons = np.linspace(-54.43, -54.40, n_horiz)
# lats = np.linspace(72.926, 72.955, n_horiz) # GLACIER C

lons = np.linspace(-54.53, -54.427, n_horiz)
lats = np.linspace(72.98, 72.967, n_horiz) # GLACIER N

# the depths
depths = np.linspace(100, 600, n_depth)

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
)


# Initial particle positions — BEFORE execute()
ax.plot(pset.lon, pset.lat, "ko",
        transform=ccrs.PlateCarree(), zorder=5, label='Start')  

#--- Run simulation ---
output_file = pset.ParticleFile(name="Upernavik_ex.zarr", outputdt=timedelta(hours=1))
pset.execute(
    parcels.AdvectionRK4,
    runtime=timedelta(days=29),
    dt=timedelta(hours=0.1),
    output_file=output_file,
)

# Final particle positions — AFTER execute()
ax.plot(pset.lon, pset.lat, "rx",
        transform=ccrs.PlateCarree(), zorder=5, label='End')    

  
def find_distance(lat2, lon2):
    
    lat1 = -54.3532829284668
    lon1 = 72.85934448242188

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


# Particle tags
tags = np.array([False for i in range(n_horiz*n_depth)])
print(tags)


# Load and plot full tracks
ds = xr.open_zarr("Upernavik_ex.zarr")
for i in range(ds.dims['trajectory']):
    
    for j in range(len(ds['lon'][i])):
        if find_distance(ds['lon'][i][j], ds['lat'][i][j])<4:
            tags[i]=True
            ax.scatter(ds['lon'][i][j], ds['lat'][i][j], transform=ccrs.PlateCarree())
            
    ax.plot(ds['lon'][i], ds['lat'][i], transform=ccrs.PlateCarree(), linewidth=1, zorder=3)


ax.plot(-54.447349548339844, 72.9875717163086, "bo",
        transform=ccrs.PlateCarree(), zorder=5)    

ax.plot(-54.36542892456055, 72.93659210205078, "bo",
        transform=ccrs.PlateCarree(), zorder=5)    

ax.plot(-54.3532829284668, 72.85934448242188, "bo",
        transform=ccrs.PlateCarree(), zorder=5, label='Outlets')  


ax.set_extent([-55.15, -54.1, 72.71, 73.05], crs=ccrs.PlateCarree())

ax.text(-55.11, 72.74, f"Depth: {depth} m", transform=ccrs.PlateCarree(), color = "white")

plt.legend(loc="lower left")
plt.show()



