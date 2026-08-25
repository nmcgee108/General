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

mat = loadmat("/Volumes/science_share/gdwbc_and_CF_sadcp_segments.mat", squeeze_me=True)
print(mat.keys())

# for k in mat.keys():
#     if not k.startswith('__'):
#         print(k, mat[k].shape, mat[k].dtype)
        
        
data = mat["vm_data"]

print(data[0, 0].dtype)

z_bin = data[0, 0]["z_bin"]

print(z_bin[0][0].squeeze().astype(float))
        

# deal with matlab times
# def matlab_datenum_to_datetime(dn):
#     return datetime.fromordinal(int(dn)) + \
#            timedelta(days=dn % 1) - \
#            timedelta(days=366)

# times = [matlab_datenum_to_datetime(t)
#          for t in mat["Xtime"].squeeze()]

print(hi)

# Helper to extract a field from the nested MATLAB struct
def get_field(mat, field):
    return mat['new2015'][field][0][0].squeeze().astype(float)

lon   = get_field(mat, 'lon')
lat   = get_field(mat, 'lat')
CT    = get_field(mat, 'CT')
SA    = get_field(mat, 'SA')
depth = get_field(mat, 'depth')

print(lon.shape, lat.shape, CT.shape, SA.shape, depth.shape)


ds = xr.Dataset(
    data_vars={
        "Conservative_Temperature": (["cast", "depth"], CT),
        "Absolute_Salinity":        (["cast", "depth"], SA),
    },
    coords={
        "depth":     depth[0, :],   # just the first cast's depth vector
        "latitude":  ("cast", lat),
        "longitude": ("cast", lon),
    }
)


ds.to_netcdf("2015_profiles.nc")
print("File Saved")










