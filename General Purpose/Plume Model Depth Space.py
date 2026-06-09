#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 09:49:48 2025

@author: nataliemcgee
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from full_plume import run_plume
from analytical_plume_gen import analytical_plume
import csv


# Load the CTD NetCDF dataset
ctd_file = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Final CTD Datasets/uc_patch_dataset_new.nc"


with Dataset(ctd_file, 'r') as ncfile:
    
    # Extract cast and depth-indexed variables (many values per cast)
    # Creates a 2D numpy array of values indexed by cast number and depth bin (missing data filled with NaN)
    depth_bin = np.ma.filled(ncfile.variables['Depth Bin'][:], np.nan)
    SA = np.ma.filled(ncfile.variables['Absolute Salinity'][:], np.nan)
    CT = np.ma.filled(ncfile.variables['Conservative Temp'][:], np.nan)

    
#-----------------------------
# THINGS TO ENTER
#-----------------------------

cast_num = 3
plume_depth = 659
Q_discharge = 60 # subglacial discharge (m3/s)
width = 500 # (m)

#-----------------------------
# Full plume model
#-----------------------------
Q0 = Q_discharge/width # (m2/s)

depth_mask = np.where(depth_bin[cast_num-1]<plume_depth)
zi = -depth_bin[cast_num-1][depth_mask][::-1] # Must input values from deepest to shallowest!
Ta = CT[cast_num-1][depth_mask][::-1]   # ambient fjord temperature at zi
Sa = SA[cast_num-1][depth_mask][::-1] # ambient fjord salinity at zi


xi = np.zeros_like(zi)
Sa1=np.average(Sa[:5])*np.ones_like(zi) # Ambient properties: average the first ?m above the plume depth
Ta1=np.average(Ta[:5])*np.ones_like(zi)
Na = np.ones_like(zi)*10 # ambient fjord nitrate at zi
alpha = 0.1 # entrainment coefficient

if max(depth_bin[cast_num-1])<plume_depth:
    print("****Warning: Chosen plume depth exceeds CTD data! T, S extrapolated at depth*****")
    zi = np.insert(zi, 0, -plume_depth)
    xi = np.insert(xi, 0, 0)
    Ta = np.insert(Ta, 0, np.average(Ta[:10]))
    Sa = np.insert(Sa, 0, np.average(Sa[:10]))
    Na = np.insert(Na, 0, 10)

sol=run_plume(zi, xi, Ta, Sa, Na, Q0, alpha)
sol1=run_plume(zi, xi, Ta1, Sa1, Na, Q0, alpha)

#-----------------------------
# Analytical plume model
#-----------------------------

S_AW=np.average(Sa[:10]) # Ambient properties: average the first ?m above the plume depth
T_AW=np.average(Ta[:10])
Tplume, Splume, AWp, SGDp, SMWp, Q_AW, Q_SMW = analytical_plume(T_AW=T_AW,
                                                                S_AW=S_AW,
                                                                Q_SGD=Q_discharge,
                                                                h_gl=plume_depth,
                                                                w=width)

#fig, axes = plt.subplots(2, 3, figsize=(12, 12))
fig, axes = plt.subplots(1, 3, figsize=(11, 7))

# axes[0,0].plot(xi, zi)
# axes[0,0].set_xlabel("Distance from Glacier Front [m]")
# axes[0,0].set_ylabel("Depth [m]")
# axes[0,0].set_title("Glacier Front Shape")

axes[0].plot(sol["b"], sol['z'], color = "darkblue", label = "Full Model (stratified)")
#axes[0].plot(sol1["b"], sol1['z'], color = "red", label = "Full Model (unstratified)")
axes[0].set_xlabel("Plume Width [m]", fontsize = 13)
axes[0].hlines(sol['zNB'], 0, max(sol["b"])+5, color = "black", linestyle = "dotted", label = "Neutral Buoyancy Depth")
axes[0].set_ylabel("Depth [m]")
#axes[0].legend()

# axes[0,2].plot(sol["u"], sol['z'], color = "darkblue")
# axes[0,2].set_xlabel("Plume Velocity [m/s]")

S_FPW = 34.25  # Fjord Polar Water (not shelf PW. See Muilwijk definitions)
T_FPW = 1.8
axes[1].plot(sol["S"], sol['z'], color = "darkblue")
axes[1].plot(sol1["S"], sol1['z'], color = "red")
axes[1].plot(sol["Sa"], sol['z'], color = "gray", label = "Ambient (stratified)", linestyle = "dotted")
axes[1].plot(sol1["Sa"], sol1['z'], color = "black", label = "Ambient (unstratified)", linestyle = "dotted")
#axes[1].vlines(S_AW, min(sol['z']), max(sol['z']), color = "gray", linestyle ="dashed", label = "Ambient (Analytic Model)")
#axes[1].scatter(sol["S"][-1], sol['z'][-1], color ="darkblue", marker = "s", label = "Full Model Output")
#axes[1].scatter(Splume, max(sol['z']), color = "darkviolet", label = "Analytic Model Output", marker = "s")
#axes[1].scatter(Splume*0.59 + S_FPW*0.41, max(sol['z']), marker = "^", s = 80, color = "darkviolet", label = "Modified Analytic Output")
axes[1].set_xlabel("Absolute Salinity [g/kg]", fontsize = 13)
axes[1].grid(True)
axes[1].set_xlim(32, 35)
# axes[1].set_ylim(-500, -250)
axes[1].legend(loc="lower left")

axes[2].plot(sol["T"], sol['z'], color = "darkblue")
axes[2].plot(sol["Ta"], sol['z'], color = "gray", linestyle = "dotted")
axes[2].plot(sol1["T"], sol1['z'], color = "red")
axes[2].plot(sol1["Ta"], sol1['z'], color = "black", linestyle = "dotted")
#axes[2].vlines(T_AW, min(sol['z']), max(sol['z']), color = "gray", linestyle ="dashed")
#axes[2].scatter(sol["T"][-1], sol['z'][-1], color ="darkblue", marker = "s")
#axes[2].scatter(Tplume, max(sol['z']), color = "darkviolet", marker = "s")
#axes[2].scatter(Tplume*0.59 + T_FPW*0.41, max(sol['z']), marker = "^", s = 80, color = "darkviolet")
axes[2].set_xlabel("Conservative Temperature [°C]", fontsize = 13)
axes[2].set_xlim(0, 3)
axes[2].grid(True)
# axes[2].set_ylim(-500, -250)


# axes[1,2].plot(sol["mdot"], sol['z'], color = "darkblue")
# axes[1,2].set_xlabel(r"Melt Rate [m$^2$/day]")

plt.suptitle(f"Plume Model Output Forced with Cast {cast_num}" + f"\nQ = {Q_discharge}"+r"m$^3$/s," + f" Width = {width}m, Plume Depth = {plume_depth}m")
plt.tight_layout()
plt.show()