#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GENERAL PURPOSE

Plot data from several CTD casts in TS space. Color code by oxygen and turbidity.

CTD data expected in NetCDF format.

Recommended environment: upernavik_env.yml

Created on Wed Jun 18 14:06:47 2025

@author: nataliemcgee
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import gsw

# Load the CTD NetCDF dataset
ctd_file = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Final CTD Datasets/uc_patch_dataset_new.nc"

with Dataset(ctd_file, 'r') as ncfile:

    #print all the variables and their dimensions:
    #print(ncfile)
    
    #print the names (keys) of the variables
    print(ncfile.variables.keys())
    
    # Extract cast-indexed variables (one value per cast)
    # Creates a 1D numpy array of values indexed by cast number (missing data filled with NaN)
    lat = np.ma.filled(ncfile.variables['Latitude'][:], np.nan) 
    lon = np.ma.filled(ncfile.variables['Longitude'][:], np.nan)
    cast_num = np.ma.filled(ncfile.variables['Cast_Number'][:], np.nan)
    date = np.ma.filled(ncfile.variables['Date'][:], np.nan)
    flag = np.ma.filled(ncfile.variables['Flag'][:], np.nan)
    
    # Extract cast and depth-indexed variables (many values per cast)
    # Creates a 2D numpy array of values indexed by cast number and depth bin (missing data filled with NaN)
    SA = np.ma.filled(ncfile.variables['Absolute Salinity'][:], np.nan)
    CT = np.ma.filled(ncfile.variables['Conservative Temp'][:], np.nan)
    turb = np.ma.filled(ncfile.variables['Turbidity'][:], np.nan)
    oxy = np.ma.filled(ncfile.variables['Oxygen'][:], np.nan)


#---------------------------------------------------------
# Set up color coding by oxygen, turbidity
#---------------------------------------------------------

colormap_oxy = plt.colormaps['viridis'] # Choose colormaps
colormap_turb = plt.colormaps['cividis']

max_oxy, max_turb = 0, 0
min_oxy, min_turb = 50, 50

# Identify (clumsily) maximum and minimum values for oxygen and turbidity
for i in range(8):
    if max(oxy[i]) > max_oxy: max_oxy = max(oxy[i])
    if max(turb[i]) > max_turb: max_turb = max(turb[i])
    if min(oxy[i]) < min_oxy: min_oxy = min(oxy[i])
    if min(turb[i]) < min_turb: min_turb = min(turb[i])
        

norm_oxy = mcolors.Normalize(vmin=min_oxy, vmax=max_oxy) # Normalize the colormaps with a max and min value
norm_turb = mcolors.Normalize(vmin=min_turb, vmax=max_turb)

sm_oxy = plt.cm.ScalarMappable(cmap=colormap_oxy, norm=norm_oxy) # Create coloring capabilities based on numerical values
sm_turb = plt.cm.ScalarMappable(cmap=colormap_turb, norm=norm_turb)

# Choose maximum and minimum values of T and S for generating isopycnal field
minS = 18
maxS = 35
minT = -2
maxT = 9

sx = np.arange(minS, maxS + 0.1, 0.1) # Create salinity array (x-axis)
ty = np.arange(minT, maxT + 0.1, 0.1) # Create temperature array (y-axis)

# Create grid of salinity and temperature
S, T = np.meshgrid(sx, ty)

# Generate a potential density field for making isopycnals
PDEN = gsw.rho(S, T, 0) - 1000  # Potential density anomaly at reference pressure = 0


#---------------------------------------------------------
# Plot TS diagrams colored by oxygen and turbidity
#---------------------------------------------------------

# Set up plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5)) 

# Draw contours (isopycnals) on both subplots
contour0 = axes[0].contour(S, T, PDEN, levels=np.arange(10, 30), colors='grey') 
contour1 = axes[1].contour(S, T, PDEN, levels=np.arange(10, 30), colors='grey')
plt.clabel(contour0, inline=True, fmt='%1.0f')  # Label contours
plt.clabel(contour1, inline=True, fmt='%1.0f')  

# Iterate through casts and plot data from each
for i in range(len(cast_num)):
    color_oxy = colormap_oxy(norm_oxy(oxy[i]))
    color_turb = colormap_turb(norm_turb(turb[i]))
    axes[0].scatter(SA[i], CT[i], c=color_oxy, s=3)
    axes[1].scatter(SA[i], CT[i], c=color_turb, s=3)

# Label the subplots
axes[0].set_title("Colored by Oxygen")
axes[0].set_xlabel("Absolute Salinity [g/kg]")
axes[0].set_ylabel("Conservative Temperature [°C]")

axes[1].set_title("Colored by Turbidity")
axes[1].set_ylabel("Conservative Temperature [°C]")
axes[1].set_xlabel("Absolute Salinity [g/kg]")

# Colorbars
cbar_oxy = fig.colorbar(sm_oxy, ax=axes[0], orientation='vertical')
cbar_oxy.set_label("Oxygen Saturation [ml/l]")

cbar_turb = fig.colorbar(sm_turb, ax=axes[1], orientation='vertical')
cbar_turb.set_label("Turbidity [NTU]")

fig.suptitle("Upernavik Fjord")
plt.tight_layout()
plt.show()


