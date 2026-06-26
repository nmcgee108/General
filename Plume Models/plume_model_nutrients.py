#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 11:31:12 2026

@author: nataliemcgee
"""

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gsw
from full_plume import run_plume
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)

# Load the CTD NetCDF dataset
ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"
nutrients_file = pd.read_csv("/Users/nataliemcgee/Documents/Upernavik Data/Nutrients/NutrientsUS2024_plotting.csv", encoding="latin-1")

ctd_ds = xr.open_dataset(ctd_netcdf)

#-----------------------------
# THINGS TO ENTER
#-----------------------------
cast_num = 3
plume_depth = 650
Q_discharge = 55 # subglacial discharge (m3/s)
width = 300 # (m)

#-----------------------------

ctd_depth = ctd_ds["depth"].values
ctd_sal = ctd_ds["SAL_ABSOLUTE"][cast_num-1].values
ctd_temp = ctd_ds["CONSERVATIVE_TEMP"][cast_num-1].values
ctd_fluor = ctd_ds["FLUORESCENCE"][cast_num-1].values
ctd_turb = ctd_ds["TURBIDITY"][cast_num-1].values
ctd_oxy = ctd_ds["OXYGEN"][cast_num-1].values
ctd_lats = ctd_ds["LAT"][cast_num-1].values
ctd_lons = ctd_ds["LON"][cast_num-1].values
ctd_castnums = ctd_ds["cast"][cast_num-1].values +1
sigma0 = gsw.sigma0(ctd_sal, ctd_temp)

valid = np.where(~np.isnan(ctd_sal))[0]

if len(valid) > 0:
    ctd_maxdepth = float(ctd_depth[valid[-1]])

# Extract nutrient data
# Using cast 8: csv lines 25-33

#cast_range = (25, 33)  # Cast 8 (downstream fjord)
#cast_range = (63, 71)   # Cast 14 (trough)
cast_range = (7, 13)    # Cast 3 (upstream fjord)

nitrate_value = pd.to_numeric(nutrients_file['NO3'][cast_range[0]:cast_range[1]])
sample_cast = nutrients_file['St#'][cast_range[0]:cast_range[1]]
sample_depth = pd.to_numeric(nutrients_file['Depth  '][cast_range[0]:cast_range[1]])

sample_depth_reversed = sample_depth.iloc[::-1]
nitrate_value_reversed = nitrate_value.iloc[::-1]

# Create a linearly interpolated nitrate profile based of cast 8 nitrate data
nitrate_profile = np.interp(ctd_depth, sample_depth_reversed, nitrate_value_reversed)
    

#-----------------------------
# Full plume model
#-----------------------------
Q0 = Q_discharge/width # (m2/s)
depth_mask = np.where(ctd_depth<plume_depth)
zi = -ctd_depth[depth_mask][::-1] # Must input values from deepest to shallowest!
xi = np.zeros_like(zi)
Ta = ctd_temp[depth_mask][::-1]   # ambient fjord temperature at zi
Sa = ctd_sal[depth_mask][::-1] # ambient fjord salinity at zi
Na = nitrate_profile[depth_mask][::-1] # ambient fjord nitrate at zi
alpha = 0.1 # entrainment coefficient

if ctd_maxdepth<plume_depth:
    print("****Warning: Chosen plume depth exceeds CTD data! T, S extrapolated at depth*****")
    zi = np.insert(zi, 0, -plume_depth)
    xi = np.insert(xi, 0, 0)
    Ta = np.insert(Ta, 0, np.average(Ta[:10]))
    Sa = np.insert(Sa, 0, np.average(Sa[:10]))
    Na = np.insert(Na, 0, np.average(Na[:10]))

sol=run_plume(zi, xi, Ta, Sa, Na, Q0, alpha)

NBD = sol["zNB"]
nitrate_NBD = sol["NNB"]
amb_nitrate_NBD = nitrate_profile[-int(sol["zNB"])]
nitrate_anomaly = nitrate_NBD-amb_nitrate_NBD
volume_flux = sol["QNB"]*width

# Units note: Multiply nitrate x 1000 to convert from liters to m^3, then
# divide by 10^6 to convert from micromoles to moles. 
NFA = (nitrate_NBD - amb_nitrate_NBD)*volume_flux / 1e3

print(f"NBD = {NBD:.2f} [meters]")
print(f"Volume Flux = {volume_flux:.2f} [m3/s]")
print(f"Nitrate @ NBD = {nitrate_NBD:.2f} [uM]")
print(f"Ambient Nitrate @ NBD = {amb_nitrate_NBD:.2f} [uM]")
print(f"Nitrate Anomaly @ NBD = {nitrate_anomaly:.2f} [uM]")
print(f"***** NFA = {NFA:.2f} [mol/s] *****")



fig, axes = plt.subplots(1,2, figsize = (10, 8))
axes[0].plot(nitrate_profile, -ctd_depth)
axes[0].set_ylim(-plume_depth, 0)
axes[0].set_xlabel(r"Shelf NO$_3$")

axes[1].plot(sol["N"], sol["z"])
axes[1].set_ylim(-plume_depth, 0)
axes[1].set_xlabel(r"Plume NO$_3$")









