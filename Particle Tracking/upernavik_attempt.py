#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 10:22:24 2026

@author: nataliemcgee
"""

import math
from datetime import timedelta

import matplotlib.pyplot as plt
import numpy as np
import trajan as ta
import xarray as xr

import cartopy.crs as ccrs
import cartopy.feature
import warnings
import parcels

warnings.filterwarnings("ignore", category=UserWarning, module="parcels")


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
print(timestamps[0])                                   # sanity check dates

fieldset = parcels.FieldSet.from_netcdf(
    filenames,
    variables,
    dimensions,
    timestamps=timestamps,
    allow_time_extrapolation=True)


# fieldset = parcels.FieldSet.from_mitgcm(filenames, 
#                                         variables, 
#                                         dimensions, 
#                                         timestamps=timestamps,
#                                         allow_time_extrapolation=True)

#print(fieldset)


# Then make sure your particle release time matches:
pset = parcels.ParticleSet.from_line(
    fieldset=fieldset,
    pclass=parcels.JITParticle,
    size=5,  # releasing 5 particles
    start=(-55.5, 73),  # releasing on a line: the start longitude and latitude
    finish=(-56.5, 73),  # releasing on a line: the end longitude and latitude
    depth = 500,
    time=start_date
)

#print(pset)

# Plot zonal velocity

fig, ax = plt.subplots(
    subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=-42)}
)

target_depth = 300 
depth_index = np.argmin(np.abs(fieldset.W.grid.depth - target_depth))
print("Selected depth:", fieldset.W.grid.depth[depth_index])

time_index = 0

fieldset.computeTimeChunk(time = 30, dt = 1)

fig, ax = plt.subplots(
    subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=-42)}

)


mesh = ax.pcolormesh(
    fieldset.W.grid.lon,
    fieldset.W.grid.lat,
    fieldset.W.data[time_index, depth_index, :, :],
    transform=ccrs.PlateCarree(),
)

fig.colorbar(mesh, ax=ax, orientation='vertical', pad=0.05)



# Plot particle positions
ax.plot(pset.lon, pset.lat, "ko")


output_file = pset.ParticleFile(
    name="Upernavik_ex.zarr", outputdt=timedelta(hours=1)
)


pset.execute(
    parcels.AdvectionRK4,
    runtime=timedelta(days=29),
    dt=timedelta(hours=4),
    output_file=output_file,
)


# Plot particle final positions
ax.plot(pset.lon, pset.lat, "rx")
print(pset)

ds = xr.open_zarr("Upernavik_ex.zarr")
ds.traj.plot(margin=2)

ax.set_extent([-55.15, -54.1, 72.71, 73.05], crs=ccrs.PlateCarree())


plt.show()











