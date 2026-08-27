#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 12:04:44 2026

@author: nataliemcgee
"""
import xarray as xr
import numpy as np
import gsw
import matplotlib.pyplot as plt
import math
import pandas as pd
import matplotlib.colors as mcolors
import numpy as np

# model_netcdf = "/Users/nataliemcgee/Downloads/NO3_202106.nc"
model_netcdf = "/Users/nataliemcgee/Desktop/Model Output/Uvel_202106.nc"
#model_netcdf = "/Users/nataliemcgee/Downloads/Theta_202106.nc"

model_ds = xr.open_dataset(model_netcdf)

# View all info
print(model_ds)
print(model_ds.depths)

depth_index = 10
depth = model_ds.depths[depth_index].values
print(depth)
    
# Pull first timestep, surface layer
for i in range(25,30):
    no3_slice = model_ds.Uvel.isel(iterations=0, depths=i).values
    lon = model_ds.longitude.values
    lat = model_ds.latitude.values
    
    plt.figure(figsize=(10, 6))
    plt.pcolormesh(lon, lat, no3_slice, cmap='viridis')
    #plt.colorbar(label='Potential Temp (°C)')
    plt.colorbar(label='Velocity')
    plt.xlim(-58, -54)
    plt.ylim(72.5, 73.5)
    depth = model_ds.depths[i].values
    plt.title(f'Along-fjord Velocity at {depth:.2f}m -- {i}')
    #plt.title(f'Temp at {depth:.2f}m -- {i}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show()