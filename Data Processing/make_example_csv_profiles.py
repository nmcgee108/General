#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 16:12:15 2026

@author: nataliemcgee
"""

import xarray as xr
import pandas as pd
import numpy as np

####### REAL PROFILE ############

ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"

ctd_ds = xr.open_dataset(ctd_netcdf)

ctd_depth = ctd_ds["depth"].values
ctd_cast = ctd_ds["cast"].values
ctd_sal = ctd_ds["SAL_ABSOLUTE"].values
ctd_temp = ctd_ds["CONSERVATIVE_TEMP"].values

data = {
    "Depth": ctd_depth[:700],
    "Absolute Salinity": ctd_sal[2][:700],
    "Temperature": ctd_temp[2][:700]
}

# Create a DataFrame
df = pd.DataFrame(data)

df.to_csv("real_profile.csv", index=False)

####### Constant/LINEAR STRAT ############

sal1 = np.linspace(32, 35, 700)
temp1 = np.linspace(-1, 3, 700)

data1 = {
    "Depth": ctd_depth[:700],
    "Absolute Salinity": sal1,
    "Temperature": temp1
}

# Create a DataFrame
df1 = pd.DataFrame(data1)

df1.to_csv("linear.csv", index=False)

####### UNSTRATIFIED ############

sal2 = np.ones(700)*34
temp2 = np.ones(700)*2

data2 = {
    "Depth": ctd_depth[:700],
    "Absolute Salinity": sal2,
    "Temperature": temp2
}

# Create a DataFrame
df2 = pd.DataFrame(data2)

df2.to_csv("unstratified.csv", index=False)




















