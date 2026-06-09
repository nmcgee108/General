#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 19:58:52 2026

Model the dependence of NBD on Q0 using the analytical model 

@author: nataliemcgee
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from analytical_plume_gen import analytical_plume
import pandas as pd
import gsw
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
    
    sigma = gsw.sigma0(Sa, Ta)

    return zi, xi, Ta, Sa, Na, sigma


zis, xis, Tas, Sas, Nas, sigmas = [], [], [], [], [], []


for file in file_list:

    df = pd.read_csv(file)

    z = df["z"]
    Ta = df["Ta"]
    Sa = df["Sa"]

    zi, xi, Ta, Sa, Na, sigma = prepare_profile(z, Ta, Sa, plume_depth)

    zis.append(zi)
    xis.append(xi)
    Tas.append(Ta)
    Sas.append(Sa)
    Nas.append(Na)
    sigmas.append(sigma)
    

z = depth_bin[cast_num-1]
Ta = CT[cast_num-1]
Sa = SA[cast_num-1]

zi, xi, Ta, Sa, Na, sigma = prepare_profile(z, Ta, Sa, plume_depth)


zis.append(zi)
xis.append(xi)
Tas.append(Ta)
Sas.append(Sa)
Nas.append(Na)
sigmas.append(sigma)



fluxes = np.arange(1, 350, 0.1)
colors = plt.cm.managua(np.linspace(0, 1, 5))
years = ["2013", "2015", "2016", "2019", "2024"]

AW_Tas = [2.51, 2.36, 2.17, 1.79, 2.7]
AW_Sas = [34.71, 34.72, 34.70, 34.69, 34.70]

PW_Tas = [0.72, 0.98, 0.26, -0.17, 1.6] #0.7]
PW_Sas = [33.98, 34.00, 34.02, 33.89, 34.13] #33.75]


plt.figure(figsize=(10, 7))

for i in range(len(years)):
    
    deep_Ta = AW_Tas[i] #np.average(Tas[i][:10])
    deep_Sa = AW_Sas[i] #np.average(Sas[i][:10])
    
    NB_depths = []
    for Q in fluxes:
        Tplume, Splume, AWp, SGDp, SMWp, Q_AW, Q_SMW = analytical_plume(T_AW=deep_Ta,
                                                                        S_AW=deep_Sa,
                                                                        Q_SGD=Q,
                                                                        h_gl=plume_depth,
                                                                        w=500)
        T_final = 0.59*Tplume +0.41*PW_Tas[i]
        S_final = 0.59*Splume +0.41*PW_Sas[i]
        sigma_final = gsw.sigma0(S_final, T_final)
        
        z = zis[i]
        sigma_prof = sigmas[i]
        
        # difference from target density
        diff = sigma_prof - sigma_final
        
        # find zero crossings
        sign_change = np.where(np.diff(np.sign(diff)))[0]
        
        if len(sign_change) > 0:
            k = sign_change[0]  # first crossing (deepest, since your array is reversed)
        
            # linear interpolation
            z1, z2 = z[k], z[k+1]
            s1, s2 = sigma_prof[k], sigma_prof[k+1]
        
            NB_depth = z1 + (sigma_final - s1) * (z2 - z1) / (s2 - s1)
        else:
            NB_depth = np.nan  # or extrapolate if you prefer

        NB_depths.append(int(NB_depth))

    NBD_smooth= savgol_filter(NB_depths, window_length=200, polyorder=5)
    
    if i == 4: 
        plt.plot(fluxes/500, NBD_smooth, color = "darkblue", linewidth =2.75, label=years[i])
    else: 
        pass
        plt.plot(fluxes/500, NBD_smooth, color = colors[i], label=years[i])
    
#plt.legend(loc="lower right", fontsize = 16)
plt.xlabel(r"Q$_0$ [m$^2$/s]", fontsize = 13)
plt.ylabel("Neutral Buoyancy Depth [m]", fontsize = 13)
plt.ylim(-400, -125)
plt.grid(True)
        
plt.show()
