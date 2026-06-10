#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 11:37:24 2026

@author: nataliemcgee
"""

import xarray as xr
import numpy as np
import gsw
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from matplotlib.colors import PowerNorm
import pandas as pd

argo_netcdf = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/argo_float_data/argo_2904018_updated.nc"

argo_ds = xr.open_dataset(argo_netcdf)

sliced_argo_ds = argo_ds.sel(TIME = slice("2024-01-01", "2026-12-31"))

temp = sliced_argo_ds["TEMP_ADJUSTED"]
pres = sliced_argo_ds["PRES_ADJUSTED"]
psal = sliced_argo_ds["PSAL_ADJUSTED"]

temp_qc = sliced_argo_ds["TEMP_ADJUSTED_QC"]
pres_qc = sliced_argo_ds["PRES_ADJUSTED_QC"]
psal_qc = sliced_argo_ds["PSAL_ADJUSTED_QC"]

SA_good = []
temp_good = []
pden_good = []

for i in range(sliced_argo_ds.sizes["TIME"]):    # iterate through each cast
    temp_cast = temp.isel(TIME=i)
    psal_cast = psal.isel(TIME=i)
    pres_cast = pres.isel(TIME=i)

    good = (temp_qc.isel(TIME=i).isin([1, 2]) &
            psal_qc.isel(TIME=i).isin([1, 2]) &
            pres_qc.isel(TIME=i).isin([1, 2])) # only plot points that are "good" or "probably good"

    lon = sliced_argo_ds["LONGITUDE"].isel(TIME=i).item()
    lat = sliced_argo_ds["LATITUDE"].isel(TIME=i).item()
    SA_cast = gsw.SA_from_SP(psal_cast, pres_cast, lon, lat)
    
    SA_cast_good = np.array(SA_cast.where(good))
    temp_cast_good = np.array(temp_cast.where(good))
    
    PDEN_cast_good = gsw.rho(SA_cast_good, temp_cast_good, 0) - 1000  # Potential density anomaly at reference pressure = 0
    PDEN_cast_good_smooth = savgol_filter(PDEN_cast_good, window_length=40, polyorder=3)

    SA_good.append(SA_cast_good)
    temp_good.append(temp_cast_good)
    pden_good.append(PDEN_cast_good_smooth)
    
    if i == 0:
        print(temp_qc)
            
section = xr.Dataset(
    data_vars={
        "Absolute Salinity": (["time", "depth"], SA_good),
        "Conservative Temperature": (["time", "depth"], temp_good),
        "Potential Density Anomaly": (["time", "depth"], pden_good)},
    coords={
        "time": sliced_argo_ds['TIME'].values,
        "depth": sliced_argo_ds['DEPTH'].values})

######################
# Plot salinity section

SA_min = float(section["Absolute Salinity"].min())
SA_max = float(section["Absolute Salinity"].max())

# Build meshgrid
X, Y = np.meshgrid(section["time"], -section["depth"])

fig, ax = plt.subplots(figsize=(12, 6))
cf = ax.contourf(X, Y, section["Absolute Salinity"].T, levels=100, norm=PowerNorm(gamma=2.5), cmap='viridis', vmin=SA_min, vmax=SA_max)
cbar = plt.colorbar(cf, ax=ax, format="%.2f")
cbar.set_label("Absolute Salinity [g/kg]", labelpad=15)

# Overlay sigma0 contours
sigma_levels = [26, 26.5, 27, 27.5, 28]
cs = ax.contour(X, Y, section["Potential Density Anomaly"].T, levels=sigma_levels, colors='white', linewidths=1)
ax.clabel(cs, inline=True, fmt="%.2f", inline_spacing = 30)

plt.tick_params(axis='both', which='major')

ax.vlines([pd.Timestamp("2024-07-31"), pd.Timestamp("2024-08-11")],
          0, -500,
          color='k',
          linestyle='dashed')

ax.set_ylabel("Depth [m]")
ax.set_xlabel("Date")

######################
# Plot temp section

CT_min = float(section["Conservative Temperature"].min())
CT_max = float(section["Conservative Temperature"].max())

# Build meshgrid
X, Y = np.meshgrid(section["time"], -section["depth"])

fig, ax = plt.subplots(figsize=(12, 6))
cf = ax.contourf(X, Y, section["Conservative Temperature"].T, levels=100, norm=PowerNorm(gamma=1), cmap='plasma', vmin=CT_min, vmax=5)
cbar = plt.colorbar(cf, ax=ax, format="%.2f")
cbar.set_label("Conservative Temperature [°C]", labelpad=15)

# Overlay sigma0 contours
sigma_levels = [26, 26.5, 27, 27.5, 28]
cs = ax.contour(X, Y, section["Potential Density Anomaly"].T, levels=sigma_levels, colors='white', linewidths=1)
ax.clabel(cs, inline=True, fmt="%.2f", inline_spacing = 30)

plt.tick_params(axis='both', which='major')


ax.vlines([pd.Timestamp("2024-07-31"), pd.Timestamp("2024-08-11")],
          0, -500,
          color='k',
          linestyle='dashed')

ax.set_ylabel("Depth [m]")
ax.set_xlabel("Date")






    