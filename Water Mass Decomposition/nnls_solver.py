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


plt.rcParams['font.size']=14

ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"

ctd_ds = xr.open_dataset(ctd_netcdf)

ctd_depth = ctd_ds["depth"].values
ctd_sal = ctd_ds["SAL_ABSOLUTE"].values
ctd_pres = ctd_ds["PRESSURE"].values
ctd_temp = ctd_ds["CONSERVATIVE_TEMP"].values
ctd_fluor = ctd_ds["FLUORESCENCE"].values
ctd_turb = ctd_ds["TURBIDITY"].values
ctd_oxy = ctd_ds["OXYGEN"].values #* 43.570 # rough convert to umol/kg
ctd_lats = ctd_ds["LAT"].values
ctd_lons = ctd_ds["LON"].values
ctd_castnums = ctd_ds["cast"].values


sigma0 = gsw.sigma0(ctd_sal, ctd_temp)
converted_oxy = ctd_oxy*44.660/((1000+sigma0)/1000)


# Enter endmember values

# Currently using fjord AW
T_AW = 2.85  # Rough estimate!
S_AW = 34.6  # Rough estimate!
O_AW = 233.10  # Rough estimate! 5.35 * 43.570 to get umol/kg

# # Currently using shelf PW
T_PW = -1.5     # Rough estimate!
S_PW = 33.2     # Rough estimate!
O_PW = 326.78   # Rough estimate! 7.5 * 43.570 to get umol/kg


# # Currently using FJORD PW
# T_PW = 0.45    # Rough estimate!
# S_PW = 33.6     # Rough estimate!
# O_PW = 270.13   # Rough estimate!


T_SMW = -87 
S_SMW = 0 
O_SMW = 1050 # From Margaret's paper!

T_SGD = 0 
S_SGD = 0 
O_SGD = 457 # From Margaret's paper!


# Build A matrix
A = np.array([
    [T_AW, T_PW, T_SMW, T_SGD],
    [S_AW, S_PW, S_SMW, S_SGD],
    [O_AW, O_PW, O_SMW, O_SGD],
    [1,    1,    1,     1    ]
])

# Statistics on data rows only (exclude last row of 1s)
Amean = np.mean(A[:-1, :], axis=1)
Astd  = np.std( A[:-1, :], ddof=0, axis=1)
Avar  = np.var( A[:-1, :2], ddof=0, axis=1)

# Normalize, then re-attach row of 1s
A_norm = (A[:-1, :] - Amean[:, np.newaxis]) / Astd[:, np.newaxis]
A_norm = np.vstack([A_norm, [1, 1, 1, 1]])

# Weights [T=121, S=97, O=1, mass=121 (same as T)]. From lindeman et al. 2024
Wx = np.array([121, 97, 10, 121])

# Apply weights and check condition number
Aw_norm = Wx[:, np.newaxis] * A_norm

print("Condition number for A:", np.linalg.cond(A))
print("Condition number for AW_Norm:", np.linalg.cond(Aw_norm))

AW_fracs = []
PW_fracs = []
SMW_fracs = []
SGD_fracs = []
residuals = []

# Pre-allocate residual arrays before the loop
RESIDUAL_T = []
RESIDUAL_S = []
RESIDUAL_O = []
RESIDUAL_M = []

cast_num = 3

for i in range(len(ctd_sal[cast_num-1])):

    # Observed values
    T_obs = ctd_temp[cast_num-1][i]
    S_obs = ctd_sal[cast_num-1][i]
    O_obs = converted_oxy[cast_num-1][i]

    d = np.array([T_obs, S_obs, O_obs, 1])

    # Normalize data point (first 3 elements), then re-attach 1
    d_norm = np.append((d[:3] - Amean) / Astd, 1)

    # Weight data (mass constraint is weighted here)
    dw_norm = Wx * d_norm

    
    if np.isnan(T_obs) or np.isnan(S_obs) or np.isnan(O_obs):
        AW_fracs.append(np.nan)
        PW_fracs.append(np.nan)
        SMW_fracs.append(np.nan)
        SGD_fracs.append(np.nan)
        residuals.append(np.nan)
        RESIDUAL_T.append(np.nan)
        RESIDUAL_S.append(np.nan)
        RESIDUAL_O.append(np.nan)
        RESIDUAL_M.append(np.nan)
        continue

    # Solve the NNLS problem
    x, residual = nnls(Aw_norm, dw_norm)


    # Residuals: observed - predicted (d - C*x)
    # RESIDUAL = dw_norm - Aw_norm @ x
    RESIDUAL = A @ x

    # Store residuals
    RESIDUAL_T.append(RESIDUAL[0])
    RESIDUAL_S.append(RESIDUAL[1])
    RESIDUAL_O.append(RESIDUAL[2])
    RESIDUAL_M.append(RESIDUAL[3])
    
    AW_fracs.append(x[0])
    PW_fracs.append(x[1])
    SMW_fracs.append(x[2])
    SGD_fracs.append(x[3])
    residuals.append(residual)

    
fig, axes = plt.subplots(2, 4, figsize=(14, 14), sharey=True) 

    
axes[0, 0].plot(ctd_temp[cast_num-1], -ctd_depth, color = "k")
axes[0, 0].plot(RESIDUAL_T, -ctd_depth, color = "gray", linewidth = 2, linestyle = "dotted")

axes[0, 1].plot(ctd_sal[cast_num-1], -ctd_depth, color = "k")
axes[0, 1].plot(RESIDUAL_S, -ctd_depth, color = "gray", linewidth = 2, linestyle = "dotted")

axes[0, 2].plot(converted_oxy[cast_num-1], -ctd_depth, color = "k")
axes[0, 2].plot(RESIDUAL_O, -ctd_depth, color = "gray", linewidth = 2, linestyle = "dotted")

axes[0, 3].plot(residuals, -ctd_depth, color = "k")

axes[0, 0].set_xlabel("CT [°C]")
axes[0, 1].set_xlabel("SA [g/kg]")
axes[0, 2].set_xlabel("Oxygen [umol/kg]")
axes[0, 3].set_xlabel("Combined Residual")

    
axes[1, 0].plot(AW_fracs, -ctd_depth, color = "C3")
axes[1, 1].plot(PW_fracs, -ctd_depth, color = "skyblue")
axes[1, 2].plot(SMW_fracs, -ctd_depth, color = "mediumpurple")
axes[1, 3].plot(SGD_fracs, -ctd_depth, color = "darkgreen")

axes[1, 0].set_xlabel("AW Fraction")
axes[1, 1].set_xlabel("PW Fraction")
axes[1, 2].set_xlabel("SMW Fraction")
axes[1, 3].set_xlabel("SGD Fraction")

axes[1, 3].set_xlim(-0.01, 0.05)


fig.suptitle(f"Cast {cast_num}")

plt.tight_layout()

plt.show()





