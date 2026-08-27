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

mat = loadmat("/Volumes/science_share/gdwbc_CF_LS_sadcp_segments.mat", squeeze_me=True)
print(mat.keys())

# for k in mat.keys():
#     if not k.startswith('__'):
#         print(k, mat[k].shape, mat[k].dtype)
        
        
data = mat["vm_data"]

print(data[0].dtype)

z_bin = data[0]["z_bin"].item().T
v_dt_bin = data[0]["v_dt_bin"].item().T
u_dt_bin = data[0]["u_dt_bin"].item().T

lat_bin = data[0]["lat_bin"].item()
lon_bin = data[0]["lon_bin"].item()

        

# deal with matlab times
# def matlab_datenum_to_datetime(dn):
#     return datetime.fromordinal(int(dn)) + \
#            timedelta(days=dn % 1) - \
#            timedelta(days=366)

# times = [matlab_datenum_to_datetime(t)
#          for t in mat["Xtime"].squeeze()]


# # Helper to extract a field from the nested MATLAB struct
# def get_field(mat, field):
#     return mat['new2015'][field][0][0].squeeze().astype(float)

# lon   = get_field(mat, 'lon')
# lat   = get_field(mat, 'lat')
# CT    = get_field(mat, 'CT')
# SA    = get_field(mat, 'SA')
# depth = get_field(mat, 'depth')

print(lon_bin.shape, lat_bin.shape, v_dt_bin.shape)


ds = xr.Dataset(
    data_vars={
        "Vvel_dt": (["location", "depth"], v_dt_bin),
        "Uvel_dt": (["location", "depth"], u_dt_bin),
    },
    coords={
        "depth": z_bin[0, :],   # just the first cast's depth vector
        "latitude":  ("location", lat_bin),
        "longitude": ("location", lon_bin),
    }
)



ds.to_netcdf("adcp_data_CF_LS.nc")
print("File Saved")










