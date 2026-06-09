#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 12:29:35 2025

@author: nataliemcgee
"""
import xarray as xr
import numpy as np
import os
from glob import glob

mar_dir = "/Users/nataliemcgee/Documents/Upernavik Project/ice_discharge/"
racmo_dir = "/Users/nataliemcgee/Documents/Upernavik Project/ice_discharge/"

mar_files = sorted(glob(os.path.join(mar_dir, "MAR*.nc")))
racmo_files = sorted(glob(os.path.join(racmo_dir, "RACMO*.nc")))

# Target stations (fjord coast IDs)
target_ids = [42704, 42733, 42773, 42968, 42972, 43247, 43649]

# Data containers
mar_data_list = []
racmo_data_list = []
date_list = []
fjord_ids, fjord_lats, fjord_lons, fjord_elevs = None, None, None, None

# Loop over MAR and RACMO files together by index
for mar_file, racmo_file in zip(mar_files, racmo_files):
    mar_ds = xr.open_dataset(mar_file)
    racmo_ds = xr.open_dataset(racmo_file)


    # Filter stations by target_ids
    station_mask = np.isin(mar_ds['station'].values, target_ids)
    if not station_mask.any():
        continue

    station_idxs = np.where(station_mask)[0]

    # Select data
    mar_dis = mar_ds["discharge"].isel(station=station_idxs)
    racmo_dis = racmo_ds["discharge"].isel(station=station_idxs)

    mar_time = mar_ds["time"]
    racmo_time = racmo_ds["time"]

    # Verify matching time dimension
    if not np.array_equal(mar_time, racmo_time):
        raise ValueError(f"Time mismatch between {os.path.basename(mar_file)} and {os.path.basename(racmo_file)}")

    # Save data
    mar_data_list.append(mar_dis)
    racmo_data_list.append(racmo_dis)
    date_list.append(mar_time)

    # Save fjord info (just once)
    if fjord_ids is None:
        fjord_ids = mar_ds["station"].isel(station=station_idxs).values
        fjord_lats = mar_ds["lat"].isel(station=station_idxs).values
        fjord_lons = mar_ds["lon"].isel(station=station_idxs).values
        fjord_elevs = mar_ds["alt"].isel(station=station_idxs).values

# Concatenate over time
dates = xr.concat(date_list, dim="time").rename({"time": "Date"})
mar_data = xr.concat(mar_data_list, dim="time").rename({"time": "Date"})
racmo_data = xr.concat(racmo_data_list, dim="time").rename({"time": "Date"})

# Wrap in DataArrays with Model dimension
mar_da = xr.DataArray(
    mar_data.T,
    dims=["Date", "Fjord"],
    coords={"Date": dates, "Fjord": fjord_ids},
    name="Discharge").expand_dims(Model=["MAR"])

racmo_da = xr.DataArray(
    racmo_data.T,
    dims=["Date", "Fjord"],
    coords={"Date": dates, "Fjord": fjord_ids},
    name="Discharge").expand_dims(Model=["RACMO"])

# Combine models
discharge_all = xr.concat([racmo_da, mar_da], dim="Model")

# Final dataset
discharge_ds = xr.Dataset(
    data_vars={"Discharge": discharge_all},
    coords={
        "Date": dates,
        "Fjord": fjord_ids,
        "Model": ["RACMO", "MAR"],
        "lat": ("Fjord", fjord_lats),
        "lon": ("Fjord", fjord_lons),
        "Depth": ("Fjord", fjord_elevs)})

# Save to NetCDF
discharge_ds.to_netcdf("/Users/nataliemcgee/Documents/Upernavik Project/ice_discharge/combined_mar_racmo_discharge1.nc")
print("Saved combined discharge dataset to 'combined_mar_racmo_discharge1.nc'")


