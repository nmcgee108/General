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
             "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Morven CTD Data 2013-2019/Morven Data - 2018.csv",
             "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Morven CTD Data 2013-2019/Morven Data - 2019.csv"]

#-----------------------------
# THINGS TO ENTER
#-----------------------------

cast_num = 3
plume_depth = 400

#-----------------------------
# Full plume model
#-----------------------------

depth_mask = np.where(depth_bin[cast_num-1]<plume_depth)
zi = -depth_bin[cast_num-1][depth_mask][::-1] # Must input values from deepest to shallowest!
Ta = CT[cast_num-1][depth_mask][::-1]   # ambient fjord temperature at zi
Sa = SA[cast_num-1][depth_mask][::-1] # ambient fjord salinity at zi

xi = np.zeros_like(zi)
Na = np.ones_like(zi)*10 # ambient fjord nitrate at zi
alpha = 0.1 # entrainment coefficient

if max(depth_bin[cast_num-1])<plume_depth:
    print("****Warning: Chosen plume depth exceeds CTD data! T, S extrapolated at depth*****")
    zi = np.insert(zi, 0, -plume_depth)
    xi = np.insert(xi, 0, 0)
    Ta = np.insert(Ta, 0, np.average(Ta[:10]))
    Sa = np.insert(Sa, 0, np.average(Sa[:10]))
    Na = np.insert(Na, 0, 10)


zis, xis, Tas, Sas, Nas = [], [], [], [], []

zis.append(zi)
xis.append(xi)
Tas.append(Ta)
Sas.append(Sa)
Nas.append(Na)


for i in range(len(file_list)):
    
    df = pd.read_csv(file_list[i])
    zi = df["z"].tolist()
    Ta = df["Ta"].tolist()
    Sa = df["Sa"].tolist()
    xi = np.zeros_like(zi)
    Na = np.ones_like(zi)*10 
    
    zis.append(zi)
    xis.append(xi)
    Tas.append(Ta)
    Sas.append(Sa)
    Nas.append(Na)

print(len(zis))

Q0s = np.arange(0, 500, 10) #Q0 is discharge volume/width. check values
colors = plt.cm.BuPu(np.linspace(0.3, 1, 5))
years = ["2013", "2015", "2016", "2018", "2019", "2024"]


plt.figure(figsize=(11, 7))

for i in range(len(years)):
    
    NB_depths = []
    for Q0 in Q0s:
        sol=run_plume(zis[i], xis[i], Tas[i], Sas[i], Na, Q0, alpha)
        NB_depths.append(sol['zNB'])
    
    plt.plot(Q0s, NB_depths, color = colors[i])
        

plt.show()