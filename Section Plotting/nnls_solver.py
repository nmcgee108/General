#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 12:03:39 2026

@author: nataliemcgee
"""

from scipy.optimize import nnls
import xarray as xr
import numpy as np
import gsw
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from matplotlib.colors import PowerNorm
import math
import pandas as pd
import matplotlib.colors as mcolors

plt.rcParams['font.size']=14

ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"

ctd_ds = xr.open_dataset(ctd_netcdf)

ctd_depth = ctd_ds["depth"].values
ctd_sal = ctd_ds["SAL_ABSOLUTE"].values
ctd_temp = ctd_ds["CONSERVATIVE_TEMP"].values
ctd_fluor = ctd_ds["FLUORESCENCE"].values
ctd_turb = ctd_ds["TURBIDITY"].values
ctd_oxy = ctd_ds["OXYGEN"].values * 43.570 # rough convert to umol/kg
ctd_lats = ctd_ds["LAT"].values
ctd_lons = ctd_ds["LON"].values
ctd_castnums = ctd_ds["cast"].values

sigma0 = gsw.sigma0(ctd_sal, ctd_temp)

# Enter endmember values

# Currently using fjord AW
T_AW = 2.75  # Rough estimate!
S_AW = 34.6  # Rough estimate!
O_AW = 233.10  # Rough estimate! 5.35 * 43.570 to get umol/kg

# Currently using shelf PW
T_PW = -1.5     # Rough estimate!
S_PW = 33.2     # Rough estimate!
O_PW = 326.78   # Rough estimate! 7.5 * 43.570 to get umol/kg

T_SMW = -87 
S_SMW = 0 
O_SMW = 1050 # From Margaret's paper!

T_SGD = 0 
S_SGD = 0 
O_SGD = 457


A = [[T_AW, T_PW, T_SMW, T_SGD],
     [S_AW, S_PW, S_SMW, S_SGD],
     [O_AW, O_PW, O_SMW, O_SGD],
     [1,    1,    1,     1]]

AW_fracs = []
PW_fracs = []
SMW_fracs = []
SGD_fracs = []

cast_num = 3

for i in range(len(ctd_sal[cast_num-1])):

    
    # Observed values
    T_obs = ctd_temp[cast_num-1][i]
    S_obs = ctd_sal[cast_num-1][i]
    O_obs = ctd_oxy[cast_num-1][i] 
    
    d = [T_obs, S_obs, O_obs, 1]
    
    if np.isnan(T_obs) or np.isnan(S_obs) or np.isnan(O_obs):
        AW_fracs.append(np.nan)
        PW_fracs.append(np.nan)
        SMW_fracs.append(np.nan)
        SGD_fracs.append(np.nan)
        continue

    # Solve the NNLS problem
    x, residual = nnls(A, d)
    
    AW_fracs.append(x[0])
    PW_fracs.append(x[1])
    SMW_fracs.append(x[2])
    SGD_fracs.append(x[3])
    
    
fig, axes = plt.subplots(2, 2, figsize=(7, 14), sharey=True) 
    
axes[0, 0].plot(AW_fracs, -ctd_depth)
axes[0, 1].plot(PW_fracs, -ctd_depth)
axes[1, 0].plot(SMW_fracs, -ctd_depth)
axes[1, 1].plot(SGD_fracs, -ctd_depth)

axes[0, 0].set_xlabel("AW Fraction")
axes[0, 1].set_xlabel("PW Fraction")
axes[1, 0].set_xlabel("SMW Fraction")
axes[1, 1].set_xlabel("SGD Fraction")




