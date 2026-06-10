#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 15:31:49 2026

@author: nataliemcgee
"""

import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math
import numpy as np
import gsw

xctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/XCTD data/converted_xctd.nc"

xctd_ds = xr.open_dataset(xctd_netcdf)

xctd_depth = xctd_ds["depth"].values
xctd_cast = xctd_ds["cast"].values
xctd_times = xctd_ds["time"].values
xctd_pract_sal = xctd_ds["Salinity"].values
xctd_cons_temp = xctd_ds["Temperature"].values
xctd_pot_temp = xctd_ds["Potential_Temperature"].values
xctd_pot_dens = xctd_ds["Potential_Density"].values
xctd_cond = xctd_ds["Conductivity"].values
xctd_lats = xctd_ds["latitude"].values
xctd_lons = xctd_ds["longitude"].values

depth2d = np.tile(xctd_depth[:, np.newaxis], (1, len(xctd_lats)))
xctd_pres = gsw.p_from_z(-depth2d, xctd_lats)
xctd_sal = gsw.SA_from_SP(xctd_pract_sal, xctd_pres, xctd_lons, xctd_lats)
print('done')



offsets = [0.00, # Cast 1
           0.03, # Cast 2
           0.04, # Cast 3
           0.02, # Cast 4
           0.07, # Cast 5
           0.00, # Cast 6
           0.05, # Cast 7
           0.02, # Cast 8
           0.03, # Cast 9
           0.03, # Cast 10
           0.00, # Cast 11
           0.04, # Cast 12
           0.02, # Cast 13
           0.00, # Cast 14
           0.01, # Cast 15
           0.00, # Cast 16
           0.03, # Cast 17
           0.02, # Cast 18
           0.02] # Cast 19

adjusted_sal = np.full_like(xctd_sal, fill_value=np.nan, dtype=float)

print(len(xctd_sal.T))

for i in range(len(xctd_sal.T)):
    adjusted_sal.T[i] = xctd_sal.T[i] + offsets[i]
    
print(adjusted_sal.T[2][45], xctd_sal.T[2][45])
    


ds = xr.Dataset(
    data_vars={
        "Absolute_Salinity": (["depth", "cast"], adjusted_sal),
        "Practical_Salinity": (["depth", "cast"], xctd_pract_sal),
        "Conservative_Temperature": (["depth", "cast"], xctd_cons_temp),
        "Potential_Temperature": (["depth", "cast"], xctd_pot_temp),
        "Potential_Density": (["depth", "cast"], xctd_pot_dens),
        "Conductivity": (["depth", "cast"], xctd_cond),
    },
    coords={
        "depth": xctd_depth,
        "cast": xctd_cast,
        "latitude": ("cast", xctd_lats),
        "longitude": ("cast", xctd_lons),
        "time": ("cast", xctd_times)
    }
)

ds.to_netcdf("adjusted_xctd.nc")

print("File Saved")










