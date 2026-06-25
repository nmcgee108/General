#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Fri Jun 13 09:23:26 2025

@author: nataliemcgee
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import math
import matplotlib.colors as mcolors
import numpy as np
import gsw

uc_patch_file = "/Users/nataliemcgee/Documents/Upernavik Data/Final CTD Datasets/uc_patch_dataset_new.nc"

ncfile = Dataset(uc_patch_file, 'r')

#print all the variables and their dimensions:
print(ncfile)

#print the names (keys) of the variables
#print(ncfile.variables.keys())

#upernavik fjord are casts 1-8, possibly 9 which is on the shelf outside the mouth

lat = ncfile.variables['Latitude'][:]
lon = ncfile.variables['Longitude'][:]
cast_num = ncfile.variables['Cast_Number'][:]
date = ncfile.variables['Date'][:]
flag = ncfile.variables['Flag'][:]


depth_bin = ncfile.variables['Depth Bin'][:]
SA = ncfile.variables['Absolute Salinity'][:]
SP = ncfile.variables['Practical Salinity'][:]
theta = ncfile.variables['Potential Temp'][:]
CT = ncfile.variables['Conservative Temp'][:]
pres = ncfile.variables['Pressure'][:]
cond = ncfile.variables['Conductivity'][:]
turb = ncfile.variables['Turbidity'][:]
fluor = ncfile.variables['Fluorescence'][:]
oxy = ncfile.variables['Oxygen'][:]


ncfile.close()


def find_distance(lat2, lon2):
    '''finds distance in km between a point and station 4 which is closest to the glacier using
    the Haversine formula'''
    
    # lat and lon of station 4:
    lat1 = lat[3]
    lon1 = lon[3]
    
    # lat and lon of a different station:
  #  lat1 = lat[9]
   # lon1 = lon[9]
        
    R = 6371.0  # Earth's radius in km

    # Convert degrees to radians
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    # Haversine formula
    a = math.sin(dphi / 2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    distance = R * c
    return distance


max_dist = find_distance(lat[7], lon[7]) #using cast 8


# Set up plot
fig, axes = plt.subplots(2, 2, figsize=(7, 14), sharey=True) #6,10

colormap = plt.colormaps['viridis']
norm = mcolors.Normalize(vmin=0, vmax=max_dist)
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)

temps = []
sals = []

for i in range(20):
    
    dist = find_distance(lat[i], lon[i])
    color = colormap(norm(dist))
    if i>=8: color = "black"

    p = gsw.p_from_z(-depth_bin[i], lat[i])
    density = gsw.sigma0(SA[i], CT[i])

    axes[0,1].plot(SA[i], -depth_bin[i], c=color)
    axes[0,0].set_ylabel("Depth [m]", fontsize = 14)
    #axes[0,0].plot(SA[i], density, c=color)
    #axes[0,0].set_ylabel("Density", fontsize = 12)
    axes[0,1].set_xlabel("Absolute Salinity [g/kg]", fontsize = 14)
    axes[0,0].axhspan(-425, -150, color = "lightskyblue", alpha = 0.06)
    axes[0,0].axhspan(-250, -130, color = "pink", alpha = 0.06)
    axes[0,1].grid(True)
    axes[0,1].set_xlim(28,35)
    axes[0,1].set_ylim(-100, 0)
    
    

    axes[0,0].plot(CT[i], -depth_bin[i], c=color)
    #axes[0,1].plot(CT[i], density, c=color)
    axes[0,0].set_xlabel("Conservative Temp [°C]", fontsize = 14)
    axes[0,1].axhspan(-425, -150, color = "lightskyblue", alpha = 0.06)
    axes[0,1].axhspan(-250, -130, color = "pink", alpha = 0.06)
    axes[0,0].grid(True)
    #axes[0,0].set_xlim(-2,4.5)
    axes[0,0].set_xlim(-2,4)


    axes[1,0].plot(turb[i], -depth_bin[i], c=color)
    #axes[1,0].plot(turb[i], density, c=color)
    axes[1,0].axhspan(-425, -150, color = "lightskyblue", alpha = 0.06)
    axes[1,0].axhspan(-250, -130, color = "pink", alpha = 0.06)
    axes[1,0].set_xlabel("Turbidity [NTU]", fontsize = 14)
    axes[1,0].grid(True)


    # axes[1,1].plot(fluor[i], -depth_bin[i], c=color)
    # axes[1,1].set_xlabel(r"Fluorescence [mg/m$^3$]")
    
    axes[1,0].set_ylabel("Depth [m]", fontsize = 14)

    #axes[1,1].plot(oxy[i], -depth_bin[i], c=color)
    axes[1,1].plot(density, -depth_bin[i], c=color)
    axes[1,1].vlines(26.98, -800, 0)
    axes[1,1].vlines(27.43, -800, 0)
    axes[1,1].axhspan(-425, -150, color = "lightskyblue", alpha = 0.06)
    axes[1,1].axhspan(-250, -130, color = "pink", alpha = 0.06)
    axes[1,1].vlines(26.95, -800, 0, color = "red")
    axes[1,1].vlines(27.2, -800, 0, color = "red")
    axes[1,1].set_xlabel("Density", fontsize = 14)
    axes[1,1].set_xlim(26.5,27.6)
    #axes[1,1].set_xlabel("Oxygen [ml/l]")
    axes[1,1].grid(True)

    
    for d, t, s in zip(depth_bin[i], CT[i], SA[i]):
        if d > 600:
            temps.append(t)
            sals.append(s)
            
print(len(turb[1]), len(CT[1]))
print("Mean Temp:", np.mean(temps), "Mean Sal:", np.mean(sals))


cbar_ax = fig.add_axes([1.03,0.1,0.05,0.8])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Distance from Innermost Station [km]', fontsize = 14)


plt.tight_layout()
plt.show()









