#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 10:32:25 2026

@author: nataliemcgee
"""

import xarray as xr
import numpy as np
from netCDF4 import Dataset


def pad_old_netcdf(old_filename, new_filename):
    
    ncfile = Dataset(old_filename, 'r')
    
    SA = ncfile.variables['Absolute Salinity']
    SP = ncfile.variables['Practical Salinity']
    CT = ncfile.variables['Conservative Temp']
    THETA = ncfile.variables['Potential Temp']
    PRES = ncfile.variables['Pressure']
    COND = ncfile.variables['Conductivity']
    OX = ncfile.variables['Oxygen']
    TURB = ncfile.variables['Turbidity']
    FLUOR = ncfile.variables['Fluorescence']
    
    depth = ncfile.variables['Depth Bin']
    
    n_casts = len(SA)
    
    max_len = max(len(SA[i])+2 for i in range(n_casts))
    
    
    
    # Create empty numpy arrays that are cast x 1309 (max depth)
        
    SA_aligned = np.full((n_casts, max_len), np.nan, dtype=np.float32)
    SP_aligned = np.full((n_casts, max_len), np.nan, dtype=np.float32)
    CT_aligned = np.full((n_casts, max_len), np.nan, dtype=np.float32)
    THETA_aligned = np.full((n_casts, max_len), np.nan, dtype=np.float32)
    PRES_aligned = np.full((n_casts, max_len), np.nan, dtype=np.float32)
    COND_aligned = np.full((n_casts, max_len), np.nan, dtype=np.float32)
    OX_aligned = np.full((n_casts, max_len), np.nan, dtype=np.float32)
    TURB_aligned = np.full((n_casts, max_len), np.nan, dtype=np.float32)
    FLUOR_aligned = np.full((n_casts, max_len), np.nan, dtype=np.float32)
    
    for i in range(n_casts):
        
        depth_i = depth[i]

        valid = depth_i >= 0
        
        depth_i = depth_i[valid]
        
        SA_i    = SA[i][valid]
        SP_i    = SP[i][valid]
        CT_i    = CT[i][valid]
        THETA_i = THETA[i][valid]
        PRES_i  = PRES[i][valid]
        COND_i  = COND[i][valid]
        OX_i    = OX[i][valid]
        TURB_i  = TURB[i][valid]
        FLUOR_i = FLUOR[i][valid]
        
        start = int(depth_i[0])         #shallowest depth
        end   = int(depth_i[-1]) + 1    #deepest depth
        
        print(start, end)
    
        SA_aligned[i, start:end] = SA_i
        SP_aligned[i, start:end] = SP_i
        CT_aligned[i, start:end] = CT_i
        THETA_aligned[i, start:end] = THETA_i
        PRES_aligned[i, start:end] = PRES_i
        COND_aligned[i, start:end] = COND_i
        OX_aligned[i, start:end] = OX_i
        TURB_aligned[i, start:end] = TURB_i
        FLUOR_aligned[i, start:end] = FLUOR_i
    
    
    ds = xr.Dataset(
        data_vars={
            "SAL_ABSOLUTE": (["cast", "depth"], SA_aligned),
            "SAL_PRACTICAL": (["cast", "depth"], SP_aligned),
            "CONSERVATIVE_TEMP": (["cast", "depth"], CT_aligned),
            "POTENTIAL_TEMP": (["cast", "depth"], THETA_aligned),
            "PRESSURE": (["cast", "depth"], PRES_aligned),
            "CONDUCTIVITY": (["cast", "depth"], COND_aligned),
            "OXYGEN": (["cast", "depth"], OX_aligned),
            "TURBIDITY": (["cast", "depth"], TURB_aligned),
            "FLUORESCENCE": (["cast", "depth"], FLUOR_aligned)
        },
        coords={
            "cast": np.arange(n_casts),
            "depth": np.arange(max_len),
        },
    )
        
    ds["LAT"] = ("cast", ncfile.variables["Latitude"][:])
    ds["LON"] = ("cast", -ncfile.variables["Longitude"][:])
    ds["DATE"] = ("cast", ncfile.variables["Date"][:])
    ds["CAST_NUM"] = ("cast", ncfile.variables["Cast_Number"][:])
    ds["FLAG"] = ("cast", ncfile.variables["Flag"][:])
    
    ncfile.close()
    
    ds.to_netcdf(new_filename)
    print(f"{new_filename} complete")


uc_patch_file = "/Users/nataliemcgee/Documents/Upernavik Data/Final CTD Datasets/uc_patch_dataset_new.nc"
dc_only_file = "/Users/nataliemcgee/Documents/Upernavik Data/Final CTD Datasets/dc_only_dataset_new.nc"
rbr_patch_file = "/Users/nataliemcgee/Documents/Upernavik Data/Final CTD Datasets/rbr_patch_dataset_new.nc"
uc_only_file = "/Users/nataliemcgee/Documents/Upernavik Data/Final CTD Datasets/uc_only_dataset_new.nc"
rbr_only_file = "/Users/nataliemcgee/Documents/Upernavik Data/Final CTD Datasets/rbr_only_dataset_new.nc"

pad_old_netcdf(uc_patch_file, "uc_patch_dataset_padded.nc")
pad_old_netcdf(dc_only_file, "dc_only_dataset_padded.nc")
pad_old_netcdf(uc_only_file, "uc_only_dataset_padded.nc")
pad_old_netcdf(rbr_only_file, "rbr_only_dataset_padded.nc")
pad_old_netcdf(rbr_patch_file, "rbr_patch_dataset_padded.nc")






    
    