#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 14:40:14 2026

@author: nataliemcgee
"""

import xarray as xr
import numpy as np
import gsw
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from matplotlib.colors import PowerNorm
import math
import pandas as pd
import matplotlib.colors as mcolors

plt.rcParams['font.size']=14

ctd_netcdf = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"
bathymetry_file = pd.read_csv("/Users/nataliemcgee/Documents/Upernavik Data/Bathymetry Data/new_fjord_bathymetry1.csv")
nutrients_file = pd.read_csv("/Users/nataliemcgee/Documents/Upernavik Data/Nutrients/NutrientsUS2024_plotting.csv", encoding="latin-1")


ctd_ds = xr.open_dataset(ctd_netcdf)

ctd_depth = ctd_ds["depth"].values
ctd_sal = ctd_ds["SAL_ABSOLUTE"].values
ctd_temp = ctd_ds["CONSERVATIVE_TEMP"].values
ctd_fluor = ctd_ds["FLUORESCENCE"].values
ctd_turb = ctd_ds["TURBIDITY"].values
ctd_oxy = ctd_ds["OXYGEN"].values
ctd_lats = ctd_ds["LAT"].values
ctd_lons = ctd_ds["LON"].values
ctd_castnums = ctd_ds["cast"].values


sigma0 = gsw.sigma0(ctd_sal, ctd_temp)


# Extract bathymetry section columns for plotting
section_distances = bathymetry_file['Distance_km']  
bed = bathymetry_file['Bed_Elevation_m'] 
section_lats = bathymetry_file['Latitude']  
section_lons = bathymetry_file['Longitude'] 

start_lon, start_lat = section_lons[0], section_lats[0]
end_lon, end_lat = section_lons.iloc[-1], section_lats.iloc[-1]

# Extract nutrient data
nitrate_value = pd.to_numeric(nutrients_file['NO3'][1:41])
sample_cast = nutrients_file['St#'][1:41]
sample_depth = -pd.to_numeric(nutrients_file['Depth  '][1:41])



def find_distance(lat1, lon1, lat2, lon2):

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


# Calculate the distance of each cast along the section and the cast max depth     
    
ctd_distances = []
ctd_maxdepths = []


for i in range(len(ctd_castnums)):   

    ctd_distances.append(find_distance(ctd_lats[i], ctd_lons[i], start_lat, start_lon))

    valid = np.where(~np.isnan(ctd_sal[i]))[0]

    if len(valid) > 0:
        ctd_maxdepths.append(float(ctd_depth[valid[-1]]))
        
    else: ctd_maxdepths.append(0)


start_cast = 2
end_cast = 9
            
section = xr.Dataset(
    data_vars={
        "Absolute Salinity": (["distance", "depth"], ctd_sal[start_cast-1:end_cast]),
        "Conservative Temperature": (["distance", "depth"], ctd_temp[start_cast-1:end_cast]),
        "Fluorescence": (["distance", "depth"], ctd_fluor[start_cast-1:end_cast]),
        "Turbidity": (["distance", "depth"], ctd_turb[start_cast-1:end_cast]),
        "Oxygen": (["distance", "depth"], ctd_oxy[start_cast-1:end_cast]),
        "Potential Density Anomaly": (["distance", "depth"], sigma0[start_cast-1:end_cast]),
        "cast_nums": (["distance"], ctd_castnums[start_cast-1:end_cast]),
        "max_depth": (["distance"], ctd_maxdepths[start_cast-1:end_cast])},
    coords={
        "distance": ctd_distances[start_cast-1:end_cast],
        "depth": ctd_depth})

section = section.sortby("distance")    # Sort by distance


########################################################################################
# Plot salinity section
########################################################################################

SA_min = 30
SA_max = 35

whole_depth_sal = False  ## ENTER WHETHER WHOLE DEPTH OR TOP SEVERAL METERS

if whole_depth_sal == True: 
    text_height = 60
    triangle_height = 40
    nitrate_max = 17.51

else: 
    text_height = 20
    triangle_height = 15
    nitrate_max = 12.5

# Build meshgrid
X, Y = np.meshgrid(section["distance"], -section["depth"])

levels = np.linspace(SA_min, SA_max, 200)

fig, ax = plt.subplots(figsize=(12, 6))
cf = ax.contourf(X, Y, section["Absolute Salinity"].T, levels=levels, norm=PowerNorm(gamma=1.3), cmap='BuPu_r', vmin=SA_min, vmax=SA_max)
cbar = plt.colorbar(cf, ax=ax, format="%.2f")
cbar.set_label("Absolute Salinity [g/kg]", labelpad=15)

# Overlay sigma0 contours
sigma_levels = [26, 26.5, 27, 27.5, 28]
smoothed_sigma = savgol_filter(section["Potential Density Anomaly"], window_length=15, polyorder=1)
cs = ax.contour(X, Y, smoothed_sigma.T, levels=sigma_levels, colors='white', linewidths=1)

label_positions = [(70, -400), (70, -100), (40, -70), (18, -30)]
ax.clabel(cs, inline=True, manual=label_positions, fmt='%.2f', inline_spacing = 20)

plt.tick_params(axis='both', which='major')

for i in range(len(section['distance'])):
    ax.vlines(section['distance'][i], -section['max_depth'][i], triangle_height, zorder = 3, color="k", linewidth = 0.5)
    ax.scatter(section['distance'][i], triangle_height, c="k", marker = "^", s = 50, zorder = 4)
    ax.text(section['distance'][i]+1, text_height, str(section['cast_nums'][i].values))
    
### Plot nutrient sample values:
    
colormap = plt.colormaps['YlGn']
norm = mcolors.Normalize(vmin=0, vmax=nitrate_max)
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)

for j in range(len(nitrate_value), 0, -1): #plot shallow samples first
    size = nitrate_value[j]*10 + 10
    color = colormap(norm(nitrate_value[j]))
    ax.scatter(ctd_distances[int(sample_cast[j])-1], sample_depth[j], color=color, edgecolor = 'k', s= size, zorder=4)

    
cbar_ax = fig.add_axes([0.95,0.1,0.02,0.78])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Nitrate [µmol N-NO3/L]', fontsize = 14)
    
#################################


ax.plot(section_distances, bed, color = "gray")

ax.fill_between(section_distances, bed, min(bed) - 50, color="gray", zorder=0)

if whole_depth_sal == True: 
    ax.text(12, -970, "Glacier")
    ax.text(98, -970, "Shelf")

    ax.set_ylim(-1000, 120)
    ax.set_xlim(0, 100)
    

else: 
    ax.set_ylim(-175, 40)
    ax.set_xlim(5, 92)


ax.invert_xaxis()
ax.set_ylabel("Depth [m]")
ax.set_xlabel("Distance Along Section [km]")

########################################################################################
# Plot temp section
########################################################################################

T_min = -2
T_max = 8

whole_depth_temp = False  ## ENTER WHETHER WHOLE DEPTH OR TOP SEVERAL METERS

if whole_depth_temp == True: 
    text_height = 60
    triangle_height = 40
    nitrate_max = 17.51

else: 
    text_height = 20
    triangle_height = 15
    nitrate_max = 12.5

        

# Build meshgrid
X, Y = np.meshgrid(section["distance"], -section["depth"])

levels = np.linspace(T_min, T_max, 300)

fig, ax = plt.subplots(figsize=(12, 6))
cf = ax.contourf(X, Y, section["Conservative Temperature"].T, levels=levels, norm=PowerNorm(gamma=0.7), cmap='RdBu_r')
cbar = plt.colorbar(cf, ax=ax, format="%.2f")
cbar.set_label("Conservative Temperature [°C]", labelpad=15)

# Overlay sigma0 contours
sigma_levels = [26, 26.5, 27, 27.5, 28]
smoothed_sigma = savgol_filter(section["Potential Density Anomaly"], window_length=15, polyorder=1)
cs = ax.contour(X, Y, smoothed_sigma.T, levels=sigma_levels, colors='dimgray', linewidths=1)

label_positions = [(70, -400), (70, -100), (40, -70), (18, -30)]
ax.clabel(cs, inline=True, manual=label_positions, fmt='%.2f', inline_spacing = 20)

plt.tick_params(axis='both', which='major')

for i in range(len(section['distance'])):
    ax.vlines(section['distance'][i], -section['max_depth'][i], triangle_height, zorder = 3, color="k", linewidth = 0.5)
    ax.scatter(section['distance'][i], triangle_height, c="k", marker = "^", s = 50, zorder = 4)
    ax.text(section['distance'][i]+1, text_height, str(section['cast_nums'][i].values))
    
    
### Plot nutrient sample values:
    
colormap = plt.colormaps['YlGn']
norm = mcolors.Normalize(vmin=0, vmax=nitrate_max)
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)

for j in range(len(nitrate_value), 0, -1): #plot shallow samples first
    size = nitrate_value[j]*10 + 10
    color = colormap(norm(nitrate_value[j]))
    ax.scatter(ctd_distances[int(sample_cast[j])-1], sample_depth[j], color=color, edgecolor = 'k', s= size, zorder=4)

    
cbar_ax = fig.add_axes([0.95,0.1,0.02,0.78])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Nitrate [µmol N-NO3/L]', fontsize = 14)
    
#################################
    
ax.plot(section_distances, bed, color = "gray")

ax.fill_between(section_distances, bed, min(bed) - 50, color="gray", zorder=0)

if whole_depth_temp == True: 
    ax.text(12, -970, "Glacier")
    ax.text(98, -970, "Shelf")

    ax.set_ylim(-1000, 120)
    ax.set_xlim(0, 100)


else: 
    ax.set_ylim(-175, 40)
    ax.set_xlim(5, 92)


ax.invert_xaxis()
ax.set_ylabel("Depth [m]")
ax.set_xlabel("Distance Along Section [km]")

########################################################################################
# Plot fluorescence section
########################################################################################

f_min = 0
f_max = 16

whole_depth_fluor = False  ## ENTER WHETHER WHOLE DEPTH OR TOP SEVERAL METERS

if whole_depth_fluor == True: 
    text_height = 60
    triangle_height = 40
    nitrate_max = 17.51

else: 
    text_height = 20
    triangle_height = 15
    nitrate_max = 12.5

        
# Build meshgrid
X, Y = np.meshgrid(section["distance"], -section["depth"])

levels = np.linspace(f_min, f_max, 300)

fig, ax = plt.subplots(figsize=(12, 6))
cf = ax.contourf(X, Y, section["Fluorescence"].T, levels=levels, norm=PowerNorm(gamma=0.6), cmap='RdYlGn')
cbar = plt.colorbar(cf, ax=ax, format="%.2f")
cbar.set_label(r"Chlorophyll Fluorescence [mg/m$^3$]", labelpad=15)

# Overlay sigma0 contours
sigma_levels = [26, 26.5, 27, 27.5, 28]
smoothed_sigma = savgol_filter(section["Potential Density Anomaly"], window_length=15, polyorder=1)
cs = ax.contour(X, Y, smoothed_sigma.T, levels=sigma_levels, colors='white', linewidths=1)

label_positions = [(70, -400), (70, -150), (40, -40), (70, -20)]
ax.clabel(cs, inline=True, manual=label_positions, fmt='%.2f', inline_spacing = 20)

plt.tick_params(axis='both', which='major')

for i in range(len(section['distance'])):
    ax.vlines(section['distance'][i], -section['max_depth'][i], triangle_height, zorder = 3, color="k", linewidth = 0.5)
    ax.scatter(section['distance'][i], triangle_height, c="k", marker = "^", s = 50, zorder = 4)
    ax.text(section['distance'][i]+1, text_height, str(section['cast_nums'][i].values))
    
    
### Plot nutrient sample values:
    
colormap = plt.colormaps['YlGn']
norm = mcolors.Normalize(vmin=0, vmax=nitrate_max)
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)

for j in range(len(nitrate_value), 0, -1): #plot shallow samples first
    size = nitrate_value[j]*10 + 10
    color = colormap(norm(nitrate_value[j]))
    ax.scatter(ctd_distances[int(sample_cast[j])-1], sample_depth[j], color=color, edgecolor = 'k', s=size, zorder=4)

    
cbar_ax = fig.add_axes([0.95,0.1,0.02,0.78])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Nitrate [µmol N-NO3/L]', fontsize = 14)
    
#################################
    
ax.plot(section_distances, bed, color = "gray")

ax.fill_between(section_distances, bed, min(bed) - 50, color="gray", zorder=0)

if whole_depth_fluor == True: 
    ax.text(12, -970, "Glacier")
    ax.text(98, -970, "Shelf")

    ax.set_ylim(-1000, 120)
    ax.set_xlim(0, 100)


else: 
    ax.set_ylim(-175, 40)
    ax.set_xlim(5, 92)


ax.invert_xaxis()
ax.set_ylabel("Depth [m]")
ax.set_xlabel("Distance Along Section [km]")

########################################################################################
# Plot turbidity section
########################################################################################

turb_min = -0.2
turb_max = 3

whole_depth_turb = False  ## ENTER WHETHER WHOLE DEPTH OR TOP SEVERAL METERS

if whole_depth_turb == True: 
    text_height = 60
    triangle_height = 40
    nitrate_max = 17.51
    turb_max = 5.5

else: 
    text_height = 20
    triangle_height = 15
    nitrate_max = 12.5
    turb_max = 3

        
# Build meshgrid
X, Y = np.meshgrid(section["distance"], -section["depth"])

levels = np.linspace(turb_min, turb_max, 300)

fig, ax = plt.subplots(figsize=(12, 6))
cf = ax.contourf(X, Y, section["Turbidity"].T, levels=levels, norm=PowerNorm(gamma=0.5), cmap='BrBG_r')
cbar = plt.colorbar(cf, ax=ax, format="%.2f")
cbar.set_label(r"Turbidity [NTU]", labelpad=15)

# Overlay sigma0 contours
sigma_levels = [26, 26.5, 27, 27.5, 28]
smoothed_sigma = savgol_filter(section["Potential Density Anomaly"], window_length=15, polyorder=1)
cs = ax.contour(X, Y, smoothed_sigma.T, levels=sigma_levels, colors='white', linewidths=1)

label_positions = [(70, -400), (70, -150), (40, -40), (70, -20)]
ax.clabel(cs, inline=True, manual=label_positions, fmt='%.2f', inline_spacing = 20)

plt.tick_params(axis='both', which='major')

for i in range(len(section['distance'])):
    ax.vlines(section['distance'][i], -section['max_depth'][i], triangle_height, zorder = 3, color="k", linewidth = 0.5)
    ax.scatter(section['distance'][i], triangle_height, c="k", marker = "^", s = 50, zorder = 4)
    ax.text(section['distance'][i]+1, text_height, str(section['cast_nums'][i].values))
    
    
### Plot nutrient sample values:
    
colormap = plt.colormaps['YlGn']
norm = mcolors.Normalize(vmin=0, vmax=nitrate_max)
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)

for j in range(len(nitrate_value), 0, -1): #plot shallow samples first
    size = nitrate_value[j]*10 + 10
    color = colormap(norm(nitrate_value[j]))
    ax.scatter(ctd_distances[int(sample_cast[j])-1], sample_depth[j], color=color, edgecolor = 'k', s=size, zorder=4)

    
cbar_ax = fig.add_axes([0.95,0.1,0.02,0.78])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Nitrate [µmol N-NO3/L]', fontsize = 14)
    
#################################
    
ax.plot(section_distances, bed, color = "gray")

ax.fill_between(section_distances, bed, min(bed) - 50, color="gray", zorder=0)

if whole_depth_turb == True: 
    ax.text(12, -970, "Glacier")
    ax.text(98, -970, "Shelf")

    ax.set_ylim(-1000, 120)
    ax.set_xlim(0, 100)


else: 
    ax.set_ylim(-175, 40)
    ax.set_xlim(5, 92)


ax.invert_xaxis()
ax.set_ylabel("Depth [m]")
ax.set_xlabel("Distance Along Section [km]")


########################################################################################
# Plot oxygen section
########################################################################################

oxy_min = 4.9
oxy_max = 12

whole_depth_ox = True  ## ENTER WHETHER WHOLE DEPTH OR TOP SEVERAL METERS

if whole_depth_ox == True: 
    text_height = 60
    triangle_height = 40
    nitrate_max = 17.51

else: 
    text_height = 20
    triangle_height = 15
    nitrate_max = 12.5
    oxy_min = 5.5

        
# Build meshgrid
X, Y = np.meshgrid(section["distance"], -section["depth"])

levels = np.linspace(oxy_min, oxy_max, 300)

fig, ax = plt.subplots(figsize=(12, 6))
cf = ax.contourf(X, Y, section["Oxygen"].T, levels=levels, norm=PowerNorm(gamma=0.6), cmap='PuOr_r')
cbar = plt.colorbar(cf, ax=ax, format="%.2f")
cbar.set_label(r"Oxygen [ml/l]", labelpad=15)

# Overlay sigma0 contours
sigma_levels = [26, 26.5, 27, 27.5, 28]
smoothed_sigma = savgol_filter(section["Potential Density Anomaly"], window_length=15, polyorder=1)
cs = ax.contour(X, Y, smoothed_sigma.T, levels=sigma_levels, colors='darkslategray', linewidths=1)

label_positions = [(70, -400), (70, -150), (40, -40), (70, -20)]
ax.clabel(cs, inline=True, manual=label_positions, fmt='%.2f', inline_spacing = 20)

plt.tick_params(axis='both', which='major')

for i in range(len(section['distance'])):
    ax.vlines(section['distance'][i], -section['max_depth'][i], triangle_height, zorder = 3, color="k", linewidth = 0.5)
    ax.scatter(section['distance'][i], triangle_height, c="k", marker = "^", s = 50, zorder = 4)
    ax.text(section['distance'][i]+1, text_height, str(section['cast_nums'][i].values))
    
    
### Plot nutrient sample values:
    
colormap = plt.colormaps['YlGn']
norm = mcolors.Normalize(vmin=0, vmax=nitrate_max)
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)

for j in range(len(nitrate_value), 0, -1): #plot shallow samples first
    size = nitrate_value[j]*10 + 10
    color = colormap(norm(nitrate_value[j]))
    ax.scatter(ctd_distances[int(sample_cast[j])-1], sample_depth[j], color=color, edgecolor = 'k', s=size, zorder=4)

    
cbar_ax = fig.add_axes([0.95,0.1,0.02,0.78])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Nitrate [µmol N-NO3/L]', fontsize = 14)
    
#################################
    
ax.plot(section_distances, bed, color = "gray")

ax.fill_between(section_distances, bed, min(bed) - 50, color="gray", zorder=0)

if whole_depth_ox == True: 
    ax.text(12, -970, "Glacier")
    ax.text(98, -970, "Shelf")

    ax.set_ylim(-1000, 120)
    ax.set_xlim(0, 100)


else: 
    ax.set_ylim(-175, 40)
    ax.set_xlim(5, 92)


ax.invert_xaxis()
ax.set_ylabel("Depth [m]")
ax.set_xlabel("Distance Along Section [km]")




    