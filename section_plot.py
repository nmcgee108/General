#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 13:36:56 2025

@author: nataliemcgee
"""

import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Geod
import numpy as np
from netCDF4 import Dataset
import gsw

# Load the CSV file into a pandas DataFrame
try:
    df = pd.read_csv('fjord_section_bathymetry.csv')
except FileNotFoundError:
    print("Error: CSV file not found.")
    exit()

# Extract columns for plotting
dist = df['Distance_km']  
bed = df['Bed_Elevation_m']  

# Add CTD lines
geod = Geod(ellps='WGS84')

# Define start and end coordinates
start_lat, start_lon = 73.055, -56.546
mid_lat, mid_lon = 72.903, -55.308
end_lat, end_lon = 72.89, -54.48

#CTD locations
latitudes = [72.91, 72.8849, 72.88805, 72.8822, 72.91523, 72.9404, 73.00466, 73.017833]
longitudes = [-55.01, -55.0707, -55.0057, -54.9242, -55.454, -55.6726, -56.08366, -56.330163]

ctd_distances = []
_, _, dist_start_to_mid = geod.inv(start_lon, start_lat, mid_lon, mid_lat)

for i in range(4):
    _, _, d = geod.inv(mid_lon, mid_lat, longitudes[i], latitudes[i]) #measure relative to center point
    ctd_distances.append((dist_start_to_mid + d)/ 1000)  # add distance from center point to start point

for i in range(4,8):
    _, _, d = geod.inv(start_lon, start_lat, longitudes[i], latitudes[i]) #measure relative to start (glacier side)
    ctd_distances.append(d / 1000)  # Convert to km
    
#distances are in order of cast number. in order of closeness it should be cast 8, 7, 6, 5, 2, (1), 3, 4
#-------------------------------------------------------
# add ctd temp data
#-------------------------------------------------------

uc_patch_file = "/Users/nataliemcgee/Documents/Upernavik Project/Upernavik Data/uc_patch_dataset_2.nc"

with Dataset(uc_patch_file, 'r') as ncfile:

    #print all the variables and their dimensions:
    #print(ncfile)
    
    lat = np.ma.filled(ncfile.variables['Latitude'][:], np.nan) 
    lon = np.ma.filled(ncfile.variables['Longitude'][:], np.nan) 
    cast_num = np.ma.filled(ncfile.variables['Cast_Number'][:], np.nan) 
    date = np.ma.filled(ncfile.variables['Date'][:], np.nan) 
    flag = np.ma.filled(ncfile.variables['Flag'][:], np.nan) 
    
    
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


first_cast = 2
last_cast = 8

# Calculate potential density to create isopycnals
sigma0 = [] 

for i in range(len(depth_bin)): #just do it for all casts we select which we want later
    sigma0.append(gsw.sigma0(SA[i], CT[i]))  # shape: (casts, depths)


#cast depths
depths = []

for i in range(last_cast):
    depths.append(-max(depth_bin[i]))
    
#conservative temperature data
ctd_temp_df = pd.DataFrame({"Depth Bin": depth_bin[first_cast-1], "Cast2": CT[first_cast-1]}) #cast 2

for i in range(first_cast, last_cast): #start at cast 3
    temp_df = pd.DataFrame({"Depth Bin": depth_bin[i], f"Cast{i+1}": CT[i]})
    ctd_temp_df = ctd_temp_df.merge(temp_df, on="Depth Bin", how="left")

    
#absolute salinity data
ctd_sa_df = pd.DataFrame({"Depth Bin": depth_bin[first_cast-1], "Cast2": SA[first_cast-1]})

for i in range(first_cast, last_cast):
    SA_df = pd.DataFrame({"Depth Bin": depth_bin[i], f"Cast{i+1}": SA[i]})
    ctd_sa_df = ctd_sa_df.merge(SA_df, on="Depth Bin", how="left")
    
    
#oxygen data
ctd_oxy_df = pd.DataFrame({"Depth Bin": depth_bin[first_cast-1], "Cast2": oxy[first_cast-1]})

for i in range(first_cast, last_cast):
    oxy_df = pd.DataFrame({"Depth Bin": depth_bin[i], f"Cast{i+1}": oxy[i]})
    ctd_oxy_df = ctd_oxy_df.merge(oxy_df, on="Depth Bin", how="left")
    
    
#turb data
ctd_turb_df = pd.DataFrame({"Depth Bin": depth_bin[first_cast-1], "Cast2": turb[first_cast-1]})

for i in range(first_cast, last_cast):
    turb_df = pd.DataFrame({"Depth Bin": depth_bin[i], f"Cast{i+1}": turb[i]})
    ctd_turb_df = ctd_turb_df.merge(turb_df, on="Depth Bin", how="left")

# Potential density anomaly data
ctd_sigma_df = pd.DataFrame({"Depth Bin": depth_bin[first_cast-1], "Cast2": sigma0[first_cast-1]})

for i in range(first_cast, last_cast):
    sig_df = pd.DataFrame({"Depth Bin": depth_bin[i], f"Cast{i+1}": sigma0[i]})
    ctd_sigma_df = ctd_sigma_df.merge(sig_df, on="Depth Bin", how="left")

#now plot the ctd data
def ctd_colormap(dataframe, ctd_distances):
    # Specify correct order of casts (matching column names)
    new_order = ["Cast8", "Cast7", "Cast6", "Cast5", "Cast2", "Cast3", "Cast4"]

    # Reorder the dataframe columns (excluding Depth Bin)
    reordered_df = dataframe[["Depth Bin"] + new_order]

    # Extract depth and salinity matrix
    depths = reordered_df["Depth Bin"].values  # Shape (M,)
    values_matrix = reordered_df.drop(columns="Depth Bin").values  # Shape (M, N)

    # Prepare coordinate grids
    distance_array = np.array([ctd_distances[i - 1] for i in [8, 7, 6, 5, 2, 3, 4]])  # Match new_order, remove cast 1
    depth_array = -np.array(depths)  # Flip to plot downward

    X, Y = np.meshgrid(distance_array, depth_array)  # Shapes (M, N)

    return X, Y, values_matrix, distance_array

cast_labels = ["8", "7", "6", "5", "2", "3", "4"]


# Plot
plt.figure(figsize=(15, 9))
X, Y, vals_matrix, distance_array = ctd_colormap(ctd_turb_df, ctd_distances)
_, _, sigma_matrix, _ = ctd_colormap(ctd_sigma_df, ctd_distances)
pc = plt.contourf(X, Y, vals_matrix, cmap='viridis', shading='gouraud')
plt.colorbar(pc, label="Turbidity [NTU]")
plt.scatter(distance_array, np.ones_like(distance_array)*15, marker = "^", color = "red")
for i in range(len(cast_labels)): plt.text(distance_array[i]-0.2, 30, cast_labels[i], fontsize=9)
plt.plot(dist, bed, color="gray")
plt.fill_between(dist, bed, min(bed)-50, color = "gray")
plt.vlines(ctd_distances, depths, 0, color = 'k')
contour = plt.contour(X, Y, sigma_matrix, levels=[26, 26.5, 27, 27.5], colors='white', linewidths=1)
plt.clabel(contour, inline=True, fontsize=9, fmt='%.1f')
plt.xlabel("Distance along section [km]") 
plt.ylabel("Bed elevation [m]")
plt.title("Turbidity Section Upernavik Fjord")
plt.text(65, -1020, "Glacier")
plt.text(1, -1020, "Mouth")
plt.tight_layout()
plt.show()

# Plot
plt.figure(figsize=(15, 9))
X, Y, vals_matrix, distance_array = ctd_colormap(ctd_sa_df, ctd_distances)
_, _, sigma_matrix, _ = ctd_colormap(ctd_sigma_df, ctd_distances)
pc = plt.contourf(X, Y, vals_matrix, cmap='viridis', shading='gouraud')
plt.colorbar(pc, label="Abolute Salinity [g/kg]")
plt.scatter(distance_array, np.ones_like(distance_array)*15, marker = "^", color = "red")
for i in range(len(cast_labels)): plt.text(distance_array[i]-0.2, 30, cast_labels[i], fontsize=9)
plt.plot(dist, bed, color="gray")
plt.fill_between(dist, bed, min(bed)-50, color = "gray")
plt.vlines(ctd_distances, depths, 0, color = 'k')
contour = plt.contour(X, Y, sigma_matrix, levels=[26, 26.5, 27, 27.5], colors='white', linewidths=1)
plt.clabel(contour, inline=True, fontsize=9, fmt='%.1f')
plt.xlabel("Distance along section [km]") 
plt.ylabel("Bed elevation [m]")
plt.title("Absolute Salinity Section Upernavik Fjord")
plt.text(65, -1020, "Glacier")
plt.text(1, -1020, "Mouth")
plt.tight_layout()
plt.show()

# Plot
plt.figure(figsize=(15, 9))
X, Y, vals_matrix, distance_array = ctd_colormap(ctd_temp_df, ctd_distances)
_, _, sigma_matrix, _ = ctd_colormap(ctd_sigma_df, ctd_distances)
pc = plt.contourf(X, Y, vals_matrix, cmap='viridis', shading='gouraud')
plt.colorbar(pc, label="Conservative Temp [°C]")
plt.scatter(distance_array, np.ones_like(distance_array)*15, marker = "^", color = "red")
for i in range(len(cast_labels)): plt.text(distance_array[i]-0.2, 30, cast_labels[i], fontsize=9)
plt.plot(dist, bed, color="gray")
plt.fill_between(dist, bed, min(bed)-50, color = "gray")
plt.vlines(ctd_distances, depths, 0, color = 'k')
contour = plt.contour(X, Y, sigma_matrix, levels=[26, 26.5, 27, 27.5], colors='white', linewidths=1)
plt.clabel(contour, inline=True, fontsize=9, fmt='%.1f')
plt.xlabel("Distance along section [km]") 
plt.ylabel("Bed elevation [m]")
plt.title("Conservative Temperature Section Upernavik Fjord")
plt.text(65, -1020, "Glacier")
plt.text(1, -1020, "Mouth")
plt.tight_layout()
plt.show()

# Plot
plt.figure(figsize=(15, 9))
X, Y, vals_matrix, distance_array = ctd_colormap(ctd_oxy_df, ctd_distances)
_, _, sigma_matrix, _ = ctd_colormap(ctd_sigma_df, ctd_distances)
pc = plt.contourf(X, Y, vals_matrix, cmap='viridis', shading='gouraud')
plt.colorbar(pc, label="Oxygen [ml/l]")
plt.scatter(distance_array, np.ones_like(distance_array)*15, marker = "^", color = "red")
for i in range(len(cast_labels)): plt.text(distance_array[i]-0.2, 30, cast_labels[i], fontsize=9)
plt.plot(dist, bed, color="gray")
plt.fill_between(dist, bed, min(bed)-50, color = "gray")
plt.vlines(ctd_distances, depths, 0, color = 'k')
contour = plt.contour(X, Y, sigma_matrix, levels=[26, 26.5, 27, 27.5], colors='white', linewidths=1)
plt.clabel(contour, inline=True, fontsize=9, fmt='%.1f')
plt.xlabel("Distance along section [km]") 
plt.ylabel("Bed elevation [m]")
plt.title("Oxygen Section Upernavik Fjord")
plt.text(65, -1020, "Glacier")
plt.text(1, -1020, "Mouth")
plt.tight_layout()
plt.show()