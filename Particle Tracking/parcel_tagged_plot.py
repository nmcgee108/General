#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 10:24:03 2026

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


#PARTICLE DEPTH
depth = 3
    
############
n_horiz = 20    # particles along the line
n_depth = 20     # particles in depth

# the line in lon/lat
lons = np.linspace(-54.43, -54.40, n_horiz)
lats = np.linspace(72.926, 72.955, n_horiz)

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
    repeatdt=timedelta(days=1)
)


# #--- Run simulation ---
# output_file = pset.ParticleFile(name="Upernavik_ex.zarr", outputdt=timedelta(hours=1))
# pset.execute(
#     parcels.AdvectionRK4,
#     runtime=timedelta(days=29),
#     dt=timedelta(hours=0.1),
#     output_file=output_file,
# )

# lon1_fixed = -54.3532829284668  # Glacier S
# lat1_fixed = 72.85934448242188

lon1_fixed = -54.447349548339844 # Glacier N
lat1_fixed = 72.9875717163086

def find_distance_vectorized(lon2, lat2, lon1=lon1_fixed, lat1=lat1_fixed):
    R = 6371.0  # km
    phi1 = np.radians(lat1)
    phi2 = np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)

    a = np.sin(dphi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(dlambda/2)**2
    c = 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))
    return R * c

ds = xr.open_zarr("Upernavik_ex.zarr")

lon = ds['lon'].values  # shape (trajectory, obs)
lat = ds['lat'].values

# Particle tags
tags = np.array([False for i in range(n_horiz*n_depth*31)])
print(tags)

dist = find_distance_vectorized(lon, lat)   # shape (trajectory, obs)
tags = np.any(dist < 5, axis=1)             # True if particle ever enters zone

print(tags)

combined_depth_grid = [depth_grid for i in range(29)]

# unique depth levels and count of True at each
unique_depths = np.unique(combined_depth_grid)
counts = np.array([np.sum(tags[combined_depth_grid == d]) for d in unique_depths])

fig, ax = plt.subplots()
ax.barh(unique_depths, (counts/n_horiz)*100)
ax.invert_yaxis()  # so shallow depths are at top, like a real depth profile
ax.set_xlabel("% of particles entering zone")
ax.set_ylabel("Depth (m)")
ax.set_title("% of C Particles that enter <1 km of N plume")
plt.show()


