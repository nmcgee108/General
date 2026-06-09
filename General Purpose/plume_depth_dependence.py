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
import pandas as pd
from scipy.signal import savgol_filter



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
Q0 = 0.12
alpha = 0.1
max_depth = 800

def prepare_profile(z, Ta, Sa, max_depth):

    z = np.asarray(z)
    Ta = np.asarray(Ta)
    Sa = np.asarray(Sa)

    # remove NaNs
    valid = ~np.isnan(Ta) & ~np.isnan(Sa)
    z, Ta, Sa = z[valid], Ta[valid], Sa[valid]

    # convert to negative depth
    zi = -z

    # truncate to plume depth
    mask = zi > -max_depth
    zi, Ta, Sa = zi[mask], Ta[mask], Sa[mask]

    # reverse so deepest → shallowest
    zi, Ta, Sa = zi[::-1], Ta[::-1], Sa[::-1]

    # extrapolate if CTD doesn't reach plume depth
    if -min(zi) < max_depth:
        zi = np.insert(zi, 0, -max_depth)
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

    zi, xi, Ta, Sa, Na = prepare_profile(z, Ta, Sa, max_depth)

    zis.append(zi)
    xis.append(xi)
    Tas.append(Ta)
    Sas.append(Sa)
    Nas.append(Na)
    

z = depth_bin[cast_num-1]
Ta = CT[cast_num-1]
Sa = SA[cast_num-1]

zi, xi, Ta, Sa, Na = prepare_profile(z, Ta, Sa, max_depth)


zis.append(zi)
xis.append(xi)
Tas.append(Ta)
Sas.append(Sa)
Nas.append(Na)


plume_depths = np.arange(-700, -300, 1)
colors = plt.cm.managua(np.linspace(0, 1, 5))
years = ["2013", "2015", "2016", "2019", "2024"]


plt.figure(figsize=(10, 7))

for i in range(len(years)):
#for i in [3,4]:
    NB_depths = []
    for plume_depth in plume_depths:
        mask = zis[i] > plume_depth
        sol=run_plume(zis[i][mask], xis[i][mask], Tas[i][mask], Sas[i][mask], Nas[i][mask], Q0, alpha)
        NB_depths.append(sol['zNB'])
    
    NBD_smooth= savgol_filter(NB_depths, window_length=20, polyorder=4)
    
    # Plot original smoothed curve
    if i == 4:
        plt.plot(plume_depths, NBD_smooth, color="darkblue", linewidth = 2.75, alpha = 0.5)
        
    else: 
        plt.plot(plume_depths, NBD_smooth, color=colors[i], alpha = 0.5)
    
    # --------- LINE OF BEST FIT ---------
    # Remove NaNs just in case
    mask = ~np.isnan(NBD_smooth)
    x_fit = plume_depths[mask][40:175]
    y_fit = NBD_smooth[mask][40:175]
    
    if len(x_fit) > 2:
        coeffs = np.polyfit(x_fit, y_fit, 1)  # linear fit
        fit_line = np.polyval(coeffs, plume_depths)
        
        if i == 4: 
            plt.plot(plume_depths[40:175], fit_line[40:175], linestyle='--', color="darkblue", linewidth = 2.75, alpha=1, zorder = 4, label=years[i])
        
        else:
            plt.plot(plume_depths[40:175], fit_line[40:175], linestyle='--', color=colors[i], alpha=1, zorder = 4, label=years[i])
        print(years[i], coeffs)
    # -----------------------------------

    
plt.legend(loc="lower right", fontsize = 16)
plt.xlabel("Plume Depth [m]", fontsize = 13)
plt.ylabel("Neutral Buoyancy Depth [m]", fontsize = 13)
plt.grid(True)
plt.ylim(-375,-225)
plt.xlim(-700,-475)
        
plt.show()
