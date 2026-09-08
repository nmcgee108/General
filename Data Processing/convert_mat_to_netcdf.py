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
print(data)

z_full = [] 
v_dt_full = [] 
u_dt_full = [] 
lat_full = []
lon_full = []
date_time = []

for i in range(len(data)):

    z_bin = data[i]["z_bin"].item().T
    v_dt_bin = data[i]["v_dt_bin"].item().T
    u_dt_bin = data[i]["u_dt_bin"].item().T
    
    lat_bin = data[i]["lat_bin"].item()
    lon_bin = data[i]["lon_bin"].item()
    
    values = data[i]["t_bin"].item()
    
    base = datetime(2026, 1, 1)  # day 1.0 = Jan 1, 00:00
    times = [base + timedelta(days=v) for v in values]
    
    print(times[0], times[-1])
    
    for z in z_bin: 
        z_full.append(z) 
    
    for v in v_dt_bin:
        v_dt_full.append(v)
        
    for u in u_dt_bin:
        u_dt_full.append(u)
        
    for lat in lat_bin:
        lat_full.append(lat)
        
    for lon in lon_bin:
        lon_full.append(lon)

    for t in times:
        date_time.append(t)
        
    

# deal with matlab times
# def matlab_datenum_to_datetime(dn):
#     return datetime.fromordinal(int(dn)) + \
#            timedelta(days=dn % 1) - \
#            timedelta(days=366)

# times = [matlab_datenum_to_datetime(t)
#          for t in mat["Xtime"].squeeze()]






ds = xr.Dataset(
    data_vars={
        "Vvel_dt": (["location", "depth"], v_dt_full),
        "Uvel_dt": (["location", "depth"], u_dt_full),
    },
    coords={
        "depth": z_bin[0, :],   # just the first cast's depth vector
        "latitude":  ("location", lat_full),
        "longitude": ("location", lon_full),
        "timestamp": ("location", date_time)
    }
)



#ds.to_netcdf("adcp_data_CF_LS.nc")
print("File Saved")










