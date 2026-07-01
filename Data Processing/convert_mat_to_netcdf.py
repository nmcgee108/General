#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 12:03:15 2026

@author: nataliemcgee
"""

        
from scipy.io import loadmat
import xarray as xr
import numpy as np
from datetime import datetime, timedelta

mat = loadmat("/Users/nataliemcgee/Downloads/MB2024_seabird_CTD.mat")
print(mat.keys())

for k in mat.keys():
    if not k.startswith('__'):
        print(k, mat[k].shape, mat[k].dtype)

#cast = np.arange(mat["cast_number"].shape[1])

# deal with matlab times
# def matlab_datenum_to_datetime(dn):
#     return datetime.fromordinal(int(dn)) + \
#            timedelta(days=dn % 1) - \
#            timedelta(days=366)

# times = [matlab_datenum_to_datetime(t)
#          for t in mat["Xtime"].squeeze()]

depth = mat["depth"][0, :].astype(np.float64) # take cast 0's depth vector, assuming all casts share it,

ds = xr.Dataset(
    data_vars={
        "Nitrate": (["cast", "depth"], mat["nitrate"]),
        # "Temperature": (["depth", "cast"], mat["TX"]),
        # "Potential_Temperature": (["depth", "cast"], mat["PTX"]),
        # "Potential_Density": (["depth", "cast"], mat["PDX"]),
        # "Conductivity": (["depth", "cast"], mat["CX"]),
    },
    coords={
        "depth": depth,   
        "cast": mat["cast_number"].squeeze(),
        "latitude": ("cast", mat["lat"].squeeze()),
        "longitude": ("cast", mat["lon"].squeeze()),
        #"time": ("cast", times)
    }
)

ds.to_netcdf("nitrate_profiles.nc")

print("File Saved")










