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

mat = loadmat("/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/XCTD data/XMB2024_binned.mat")
print(mat.keys())

for k in mat.keys():
    if not k.startswith('__'):
        print(k, mat[k].shape, mat[k].dtype)
        

depth = mat["ZX"].squeeze()
cast = np.arange(mat["SX"].shape[1])


# deal with matlab times
def matlab_datenum_to_datetime(dn):
    return datetime.fromordinal(int(dn)) + \
           timedelta(days=dn % 1) - \
           timedelta(days=366)

times = [matlab_datenum_to_datetime(t)
         for t in mat["Xtime"].squeeze()]

ds = xr.Dataset(
    data_vars={
        "Salinity": (["depth", "cast"], mat["SX"]),
        "Temperature": (["depth", "cast"], mat["TX"]),
        "Potential_Temperature": (["depth", "cast"], mat["PTX"]),
        "Potential_Density": (["depth", "cast"], mat["PDX"]),
        "Conductivity": (["depth", "cast"], mat["CX"]),
    },
    coords={
        "depth": depth,
        "cast": cast,
        "latitude": ("cast", mat["Xlat"].squeeze()),
        "longitude": ("cast", mat["Xlon"].squeeze()),
        "time": ("cast", times)
    }
)

ds.to_netcdf("converted_profiles.nc")

print("File Saved")










