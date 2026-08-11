#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 08:54:35 2026

@author: nataliemcgee
"""

import xarray as xr
import numpy as np
import gsw
import matplotlib.pyplot as plt

plt.rcParams['font.size']=18

nutrient_profiles = "/Users/nataliemcgee/Documents/Upernavik Data/Nutrients/nitrate_profiles.nc"
ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"


nitrate_ds = xr.open_dataset(nutrient_profiles)
ctd_ds = xr.open_dataset(ctd_netcdf)

ctd_sal = ctd_ds["SAL_ABSOLUTE"].values
ctd_temp = ctd_ds["CONSERVATIVE_TEMP"].values
ctd_castnums = ctd_ds["cast"].values

sigma0 = gsw.sigma0(ctd_sal, ctd_temp)

depth = nitrate_ds["depth"].values
nitrate = nitrate_ds["Nitrate"].values
castnums = nitrate_ds["cast"].values


# Pad the ctd data so all casts are the same length as nitrate cast data
target_len = nitrate.shape[1]

sal_padded = np.pad(ctd_sal, ((0, 0), (0, target_len - ctd_sal.shape[1])), constant_values=np.nan)
temp_padded = np.pad(ctd_temp, ((0, 0), (0, target_len - ctd_sal.shape[1])), constant_values=np.nan)
sigma0_padded = np.pad(sigma0, ((0, 0), (0, target_len - sigma0.shape[1])), constant_values=np.nan)

sigma0_scale = np.linspace(12.3, 27.7, target_len) # create an evenly spaced density range to interpolate onto

fig, axes = plt.subplots(1, 1, figsize=(8, 10))

mask = ~np.isnan(sigma0_padded[14]) & ~np.isnan(nitrate[14])
sorted_sigma0_shelf, sorted_nitrate_shelf = zip(*sorted(zip(sigma0_padded[14][mask], nitrate[14][mask])))
nitrate_sigma0_interp_shelf = np.interp(sigma0_scale, sorted_sigma0_shelf, sorted_nitrate_shelf)


for i in [1, 2]:
    
    if np.isin(castnums[i], [1, 5, 10, 11, 12]):
        continue
    
    print(castnums[i])
    zorder = 1
    linewidth = 2
    label = None
    
    mask = ~np.isnan(sigma0_padded[i]) & ~np.isnan(nitrate[i])
    sorted_sigma0, sorted_nitrate = zip(*sorted(zip(sigma0_padded[i][mask], nitrate[i][mask])))
    nitrate_sigma0_interp = np.interp(sigma0_scale, sorted_sigma0, sorted_nitrate)

    axes.plot(nitrate_sigma0_interp - nitrate_sigma0_interp_shelf, sigma0_scale)

axes.legend(loc="lower left")
axes.set_xlabel(r"Nitrate Anomaly [$\mu$M]")


# DENSITY PROFILE STUFF
axes.invert_yaxis()
axes.set_ylabel(r"Potential Density Anomaly [kg/m$^3$]")
axes.vlines(0, 27.8, 26.2, color="k")
axes.set_ylim(27.7, 26.25)
axes.set_xlim(-5, 12)
axes.axhspan(27.05, 27.45, color = 'skyblue', alpha = 0.3)
axes.axhspan(25.7, 26.8, color = 'thistle', alpha = 0.3)







