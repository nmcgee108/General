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

mat = loadmat("/Users/nataliemcgee/Downloads/MB2026CTDdowncasts.mat", squeeze_me=True)

print(mat.keys())

# for k in mat.keys():
#     if not k.startswith('__'):
#         print(k, mat[k].shape, mat[k].dtype)
        

# deal with matlab times
# def matlab_datenum_to_datetime(dn):
#     return datetime.fromordinal(int(dn)) + \
#            timedelta(days=dn % 1) - \
#            timedelta(days=366)

# times = [matlab_datenum_to_datetime(t)
#          for t in mat["Xtime"].squeeze()]

fluor = mat["fluor"]
oxy = mat["oxy"]
pres = mat["pres"]
sal = mat["sal"]
temp = mat["temp"]
turb = mat["turb"]

lat = mat["lat"]
lon = mat["lon"]
date = mat["date"]
time = mat["utctime"]

depth = mat["depth"]

ds = xr.Dataset(
    data_vars={
        "Fluor": (["location", "depth"], fluor),
        "Oxy": (["location", "depth"], oxy),
        "Pres": (["location", "depth"], pres),
        "Sal": (["location", "depth"], sal),
        "Temp": (["location", "depth"], temp),
        "Turb": (["location", "depth"], turb)
    },
    coords={
        "depth": depth[0, :],   # just the first cast's depth vector
        "latitude":  ("location", lat),
        "longitude": ("location", lon),
        "date": ("location", date)
    }
)



ds.to_netcdf("2026_upernavik_ctd.nc")
print("File Saved")










