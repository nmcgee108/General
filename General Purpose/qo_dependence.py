#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 17:21:51 2026

@author: nataliemcgee
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from full_plume import run_plume
from scipy.signal import savgol_filter
import pandas as pd



# Load the CTD NetCDF dataset
ctd_file = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Final CTD Datasets/uc_patch_dataset_new.nc"


with Dataset(ctd_file, 'r') as ncfile:
    
    # Extract cast and depth-indexed variables (many values per cast)
    # Creates a 2D numpy array of values indexed by cast number and depth bin (missing data filled with NaN)
    depth_bin = np.ma.filled(ncfile.variables['Depth Bin'][:], np.nan)
    SA = np.ma.filled(ncfile.variables['Absolute Salinity'][:], np.nan)
    CT = np.ma.filled(ncfile.variables['Conservative Temp'][:], np.nan)
    
file_list = [ "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Morven CTD Data 2013-2019/Morven Data - 2013.csv",
             "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Morven CTD Data 2013-2019/Morven Data - 2015.csv",
             "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Morven CTD Data 2013-2019/Morven Data - 2016.csv",
             "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Morven CTD Data 2013-2019/Morven Data - 2019.csv"]

#-----------------------------
# THINGS TO ENTER
#-----------------------------

cast_num = 3
plume_depth = 650
alpha = 0.1


def prepare_profile(z, Ta, Sa, plume_depth):

    z = np.asarray(z)
    Ta = np.asarray(Ta)
    Sa = np.asarray(Sa)

    # remove NaNs
    valid = ~np.isnan(Ta) & ~np.isnan(Sa)
    z, Ta, Sa = z[valid], Ta[valid], Sa[valid]

    # convert to negative depth
    zi = -z

    # truncate to plume depth
    mask = zi > -plume_depth
    zi, Ta, Sa = zi[mask], Ta[mask], Sa[mask]

    # reverse so deepest → shallowest
    zi, Ta, Sa = zi[::-1], Ta[::-1], Sa[::-1]

    # extrapolate if CTD doesn't reach plume depth
    if -min(zi) < plume_depth:
        print("****Warning: plume depth exceeds CTD data*****")

        zi = np.insert(zi, 0, -plume_depth)
        Ta = np.insert(Ta, 0, np.mean(Ta[:10]))
        Sa = np.insert(Sa, 0, np.mean(Sa[:10]))

    xi = np.zeros_like(zi)
    Na = np.ones_like(zi) * 10

    return zi, xi, Ta, Sa, Na


zis, xis, Tas, Sas, Nas = [], [], [], [], []


for file in file_list:

    df = pd.read_csv(file)

    z = df["z"]
    Ta = df["Ta"]
    Sa = df["Sa"]

    zi, xi, Ta, Sa, Na = prepare_profile(z, Ta, Sa, plume_depth)

    zis.append(zi)
    xis.append(xi)
    Tas.append(Ta)
    Sas.append(Sa)
    Nas.append(Na)
    

z = depth_bin[cast_num-1]
Ta = CT[cast_num-1]
Sa = SA[cast_num-1]

zi, xi, Ta, Sa, Na = prepare_profile(z, Ta, Sa, plume_depth)


zis.append(zi)
xis.append(xi)
Tas.append(Ta)
Sas.append(Sa)
Nas.append(Na)


fluxes = np.arange(1, 350, 1)
Q0s = fluxes/500 #Q0 is discharge volume/width. check values
colors = plt.cm.managua(np.linspace(0, 1, 5))
years = ["2013", "2015", "2016", "2019", "2024"]


plt.figure(figsize=(10, 7))

for i in range(len(years)-1):
    
    NB_depths = []
    for Q0 in Q0s:
        sol=run_plume(zis[i], xis[i], Tas[i], Sas[i], Nas[i], Q0, alpha)
        NB_depths.append(sol['zNB'])
        
    NBD_smooth= savgol_filter(NB_depths, window_length=20, polyorder=4)
    
    #plt.plot(Q0s, NBD_smooth, color = colors[i], label=years[i])
    
#plotting 2024 seperately
    
NB_depths = []
for Q0 in Q0s:
    sol=run_plume(zis[-1], xis[-1], Tas[-1], Sas[-1], Nas[-1], Q0, alpha)
    NB_depths.append(sol['zNB'])
    
NBD_smooth= savgol_filter(NB_depths, window_length=20, polyorder=4)

plt.plot(Q0s, NBD_smooth, color = "darkblue", linewidth = 2.75, label=years[-1])
    
#plt.legend(loc="lower right", fontsize = 16)
plt.xlabel(r"Q$_0$ [m$^2$/s]", fontsize = 13)
plt.ylabel("Neutral Buoyancy Depth [m]", fontsize = 13)
plt.grid(True)
plt.ylim(-400, -125)
plt.show()
