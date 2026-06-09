#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GENERAL PURPOSE

Plot data from several CTD casts in depth space showing a variety of variables. 
Color code by distance from a point.

CTD data expected in NetCDF format.

Recommended environment: upernavik_env.yml

Created on Mon Jun 16 17:16:03 2025

@author: nataliemcgee
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pyproj import Geod
import numpy as np

# Load the CTD NetCDF dataset
ctd_file = "/Users/nataliemcgee/Documents/Upernavik Project/Upernavik Data/uc_patch_dataset_2.nc"

with Dataset(ctd_file, 'r') as ncfile:

    # Print all the variables and their dimensions:
    #print(ncfile)
    
    # Print the names (keys) of the variables
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
    depth_bin = np.ma.filled(ncfile.variables['Depth Bin'][:], np.nan)
    SA = np.ma.filled(ncfile.variables['Absolute Salinity'][:], np.nan)
    SP = np.ma.filled(ncfile.variables['Practical Salinity'][:], np.nan)
    theta = np.ma.filled(ncfile.variables['Potential Temp'][:], np.nan)
    CT = np.ma.filled(ncfile.variables['Conservative Temp'][:], np.nan)
    pres = np.ma.filled(ncfile.variables['Pressure'][:], np.nan)
    cond = np.ma.filled(ncfile.variables['Conductivity'][:], np.nan)
    turb = np.ma.filled(ncfile.variables['Turbidity'][:], np.nan)
    fluor = np.ma.filled(ncfile.variables['Fluorescence'][:], np.nan)
    oxy = np.ma.filled(ncfile.variables['Oxygen'][:], np.nan)



#---------------------------------------------------------
# Set up color coding by distance
#---------------------------------------------------------


# Define the geod and desired ellipsoid. WGS84 is standard for Earth
# Used to calculate great circle distances between lat, lon points 
geod = Geod(ellps='WGS84')

lon_ref, lat_ref = 55.01 , 72.91
distances = [] # Distances will fill the list indexed by cast order

# Calculate distance between every point and reference point
for i in range(len(cast_num)):
    _, _, d = geod.inv(lon_ref, lat_ref, lon[i], lat[i]) # geod.inv returns three arguments; only the third (distance, m) is useful 
    distances.append(d / 1000)  # Convert to km


max_dist = max(distances) 
colormap = plt.colormaps['viridis'] # Choose colormap
norm = mcolors.Normalize(vmin=0, vmax=max_dist)  # Normalize the colormap with a max and min value
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)   # Creates coloring capabilities based on numerical values


#---------------------------------------------------------
# Plot color coded cast data
#---------------------------------------------------------

# Set up plot
fig, axes = plt.subplots(2, 4, figsize=(10, 10), sharey=True, constrained_layout=True) #Constrained_layout helps format colorbar

# Iterate through casts and plot each variable in a different plot
for i in range(len(cast_num)):
    
    dist = distances[i]
    color = colormap(norm(dist)) # Assigns each cast a color in the spectrum
    
    axes[0,0].plot(SA[i], -depth_bin[i], c=color)
    axes[0,0].set_ylabel("Depth [m]")
    axes[0,0].set_xlabel("Absolute Salinity [g/kg]")

    axes[0,1].plot(SP[i], -depth_bin[i], c=color)
    axes[0,1].set_xlabel("Practical Salinity [PSU]")

    axes[0,2].plot(theta[i], -depth_bin[i], c=color)
    axes[0,2].set_xlabel("Potential Temp [°C]")

    axes[0,3].plot(CT[i], -depth_bin[i], c=color)
    axes[0,3].set_xlabel("Conservative Temp [°C]")

    axes[1,0].plot(cond[i], -depth_bin[i], c=color)
    axes[1,0].set_xlabel("Conductivity [mS/cm]")
    axes[1,0].set_ylabel("Depth [m]")

    axes[1,1].plot(turb[i], -depth_bin[i], c=color)
    axes[1,1].set_xlabel("Turbidity [NTU]")

    axes[1,2].plot(fluor[i], -depth_bin[i], c=color)
    axes[1,2].set_xlabel(r"Fluorescence [mg/m$^3$]")

    axes[1,3].plot(oxy[i], -depth_bin[i], c=color)
    axes[1,3].set_xlabel("Oxygen [ml/l]")

fig.suptitle("Upernavik Fjord")
cbar_ax = fig.add_axes([1.03,0.1,0.05,0.8])
cbar = fig.colorbar(sm, cax=cbar_ax)  # Add a colorbar
cbar.set_label('Distance from Innermost Station [km]') # Label colorbar


plt.show()







