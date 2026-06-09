#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 15:57:57 2025

@author: nataliemcgee
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pyproj import Geod
import numpy as np
import gsw
import matplotlib.patches as mpatches
from matplotlib import colors

plt.rcParams['font.size'] = 14

uc_patch_file = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Final CTD Datasets/uc_patch_dataset_new.nc"

ncfile = Dataset(uc_patch_file, 'r')

#print all the variables and their dimensions:
#print(ncfile)

#print the names (keys) of the variables
#print(ncfile.variables.keys())

#upernavik fjord are casts 1-8, possibly 9 which is on the shelf outside the mouth

#NOTE: indexes are as such: CT[7] is cast 8

lat = ncfile.variables['Latitude'][:]
lon = ncfile.variables['Longitude'][:]
cast_num = ncfile.variables['Cast_Number'][:]
date = ncfile.variables['Date'][:]
flag = ncfile.variables['Flag'][:]


SA = ncfile.variables['Absolute Salinity'][:]
CT = ncfile.variables['Conservative Temp'][:]
turb = ncfile.variables['Turbidity'][:]
oxy = ncfile.variables['Oxygen'][:]
depth = ncfile.variables['Depth Bin'][:]


ncfile.close()


colormap_oxy = plt.colormaps['viridis']
colormap_turb = plt.colormaps['cividis']

max_oxy, max_turb = 0, 0

min_oxy, min_turb = 50, 50


for i in range(8):
    if max(oxy[i]) > max_oxy: max_oxy = max(oxy[i])
    if max(turb[i]) > max_turb: max_turb = max(turb[i])
    if min(oxy[i]) < min_oxy: min_oxy = min(oxy[i])
    if min(turb[i]) < min_turb: min_turb = min(turb[i])
        

norm_oxy = mcolors.Normalize(vmin=min_oxy, vmax=max_oxy)
norm_turb = mcolors.Normalize(vmin=min_turb, vmax=max_turb)
#norm_depth = mcolors.Normalize(vmin=-800, vmax=0)
norm_depth = colors.PowerNorm(gamma=1, vmin=-800, vmax=0)
norm_castnum = mcolors.Normalize(vmin=10, vmax=19)

sm_oxy = plt.cm.ScalarMappable(cmap=colormap_oxy, norm=norm_oxy)
sm_turb = plt.cm.ScalarMappable(cmap=colormap_turb, norm=norm_turb)
sm_depth = plt.cm.ScalarMappable(cmap=colormap_oxy, norm=norm_depth)
sm_castnum = plt.cm.ScalarMappable(cmap=colormap_oxy, norm=norm_castnum)

fig, axes = plt.subplots(1, 2, figsize=(19, 8))

minS = 32
maxS = 35
minT = -2
maxT = 4

sx = np.arange(minS, maxS + 0.1, 0.1)
ty = np.arange(minT, maxT + 0.1, 0.1)


# Create grid of salinity and temperature
S, T = np.meshgrid(sx, ty)
PDEN = gsw.rho(S, T, 0) - 1000  # Potential density anomaly at reference pressure = 0

#mixing line params
T_deepwater = 2.550228
S_deepwater = 34.697933
T_ice_eff = -90
    
def runoff_line(x):
    T_deepwater = 2.550228
    S_deepwater = 34.697933
    slope = T_deepwater/S_deepwater
    return slope*x

def melting_line(x):
    slope = (T_deepwater-T_ice_eff)/S_deepwater
    return slope*x+T_ice_eff

########### PLOT 1 ########
# Plot the filled contour between 26.5 and 27.0
filled = axes[0].contourf(
    S, T, PDEN,
    levels=[26.9, 27.4],  # Just the range you want shaded
    colors=['tan'],  # Or use a colormap like cmap='Greys'
    alpha=0.5,             # Optional: make it slightly transparent
    extend='neither',       # Controls if arrows show on the ends
    )

# Plot the isopycnal contour lines
contour = axes[0].contour(S, T, PDEN, levels=[25.5, 26, 26.5, 27, 27.5, 28], colors='grey')
contour1= axes[0].contour(S, T, PDEN, levels=[27.2], linestyles = "dashed", colors='brown')

# Manually label contours
label_positions = [(32.8, 3.1), (33.1, 2.5), (33.8, 2.9), (34.5, 1.4), (35.0, 1)]
plt.clabel(contour, inline=True, manual=label_positions, fmt='%1.1f')


bottom_left_x = 2
bottom_left_y = 3
width = 0.05
height = 0.81

# # Create the rectangle patch
# rectangle = mpatches.Rectangle((34.69, 1.7), width, height,
#                                facecolor='lightblue', edgecolor='blue', linewidth=2, alpha=0.7, label = "AW (Mortensen 2021)")

# # Add the rectangle to the axes
# axes[0].add_patch(rectangle)

width_pw = 0.13
height_pw = 2.15

# # Create the rectangle patch
# rectangle_pw = mpatches.Rectangle((33.89, -0.17), width_pw, height_pw,
#                                facecolor='thistle', edgecolor='purple', linewidth=2, alpha=0.7, label = "PW (Mortensen 2021)")

# # Add the rectangle to the axes
# axes[0].add_patch(rectangle_pw)

t400 = False
t300 = False
t500 = False
t150 = False
t125 = False

#for i in range(8,19):   
for i in range(8,19):  
    color = colormap_oxy(norm_castnum(i))
    #color = colormap_oxy(norm_depth(-depth[i]))
    color = 'darkgray'
    axes[0].scatter(SA[i], CT[i], c=color, s = 3)
    
    for j in range(len(depth[i])):
        if float(depth[i][j])==500: 
            axes[0].scatter(SA[i][j], CT[i][j], color ="deeppink", s = 10, zorder = 3)
            if t500 == False: 
                axes[0].text(SA[i][j]+0.05, CT[i][j]-0.2, "500m", color = "deeppink")
                t500=True
    
    # for j in range(len(depth[i])):
    #     if float(depth[i][j])==400: 
    #         axes[0].scatter(SA[i][j], CT[i][j], color ="darkmagenta", s = 10, zorder = 3)
    #         if t400 == False: 
    #             axes[0].text(SA[i][j]+0.03, CT[i][j]-0.1, "400m", color = "darkmagenta")
    #             t400=True
                
    for j in range(len(depth[i])):
        if float(depth[i][j])==300: 
            axes[0].scatter(SA[i][j], CT[i][j], color ="blueviolet", s = 10, zorder = 3)
            if t300 == False: 
                axes[0].text(SA[i][j]-0.1, CT[i][j]+0.8, "300m", color = "blueviolet")
                t300=True
    
    for j in range(len(depth[i])):
        if float(depth[i][j])==110: 
            axes[0].scatter(SA[i][j], CT[i][j], color ="blue", s = 10, zorder = 3)
            if t150 == False: 
                axes[0].text(SA[i][j]-0.2, CT[i][j]+0.4, "150m", color = "blue")
                t150=True
                
    # for j in range(len(depth[i])):
    #     if float(depth[i][j])==125: 
    #         axes[0].scatter(SA[i][j], CT[i][j], color ="cyan", s = 10, zorder = 3)
    #         if t125 == False: 
    #             axes[0].text(SA[i][j]-0.2, CT[i][j]+0.4, "125m", color = "cyan")
    #             t125=True
                
#axes[0].plot(SA[8], CT[8], c='k', label="Cast 9 (near fjord)")
# axes[0].plot(SA[32], CT[32], c='red')
# axes[0].plot(SA[31], CT[31], c='red')
# axes[0].plot(SA[30], CT[30], c='red', label = "Cast 31 (far shelf)")
# axes[0].plot([33.25, 34.65], [-1.5, 3.25], color = "red", linestyle="dashdot")

#axes[0].set_title(r"$\mathbf{Trough/Shelf:}$" + "\nColored by Cast Number")
axes[0].set_title(r"$\mathbf{Trough/Shelf:}$")
axes[0].set_xlabel("Absolute Salinity [g/kg]")
axes[0].set_ylabel("Conservative Temperature [°C]")
axes[0].set_xlim(33.0, 35)
axes[0].set_ylim(-2, 4)


axes[0].plot(sx, runoff_line(sx), color ="orange", linestyle = "dashed", 
             #label = "Runoff Line", 
             zorder = 0)
axes[0].plot(sx, melting_line(sx), color ="pink", linestyle = "dashed", 
             #label = "Melting Line", 
             zorder = 0)

axes[0].legend()

# Colorbars

#cbar_depth = fig.colorbar(sm_castnum, ax=axes[0], orientation='vertical')
#cbar_depth = fig.colorbar(sm_depth, ax=axes[0], orientation='vertical')
#cbar_depth.set_label("Cast Number")


#---------------------------------------------------------
# Set up color coding by distance
#---------------------------------------------------------


# Define the geod and desired ellipsoid. WGS84 is standard for Earth
# Used to calculate great circle distances between lat, lon points 
geod = Geod(ellps='WGS84')

#lon_ref, lat_ref = 54.4461, 72.9971
lon_ref, lat_ref = lon[9],lat[9] #alternative distance reference
distances = [] # Distances will fill the list indexed by cast order

# Calculate distance between every point and reference point
for i in range(25):
    _, _, d = geod.inv(lon_ref, lat_ref, lon[i], lat[i]) # geod.inv returns three arguments; only the third (distance, m) is useful 
    distances.append(d / 1000)  # Convert to km


max_dist = max(distances) 
min_dist = min(distances) 
colormap = plt.colormaps['viridis'] # Choose colormap
norm = mcolors.Normalize(vmin=min_dist, vmax=max_dist)  # Normalize the colormap with a max and min value
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)   # Creates coloring capabilities based on numerical values


######### PLOT 2 ########
# Plot the filled contour between 26.5 and 27.0
filled = axes[1].contourf(
    S, T, PDEN,
    levels=[26.98, 27.43],  # Just the range you want shaded
    colors=['lightskyblue'],  # Or use a colormap like cmap='Greys'
    alpha=0.3,             # Optional: make it slightly transparent
    extend='neither'       # Controls if arrows show on the ends
    )

# filled = axes[1].contourf(
#     S, T, PDEN,
#     levels=[26.95, 27.2],
#     colors=['pink'],  # Or use a colormap like cmap='Greys'
#     alpha=0.3,             # Optional: make it slightly transparent
#     extend='neither'       # Controls if arrows show on the ends
#     )


# Plot the isopycnal contour lines
contour = axes[1].contour(S, T, PDEN, levels=[25.5, 26, 26.5, 27, 27.5, 28], colors='grey')
#contour1= axes[1].contour(S, T, PDEN, levels=[27.2], linestyles = "dashed", colors='brown')

# Manually label contours
label_positions = [(32.8, 3.1), (33.1, 2.5), (33.8, 2.9), (34.5, 1.4), (35.0, 1)]
plt.clabel(contour, inline=True, manual=label_positions, fmt='%1.1f')



width = 0.05
height = 0.81

# Create the rectangle patch
# rectangle = mpatches.Rectangle((34.69, 1.7), width, height,
#                                facecolor='lightblue', edgecolor='blue', linewidth=2, alpha=0.7)

# Add the rectangle to the axes
#axes[1].add_patch(rectangle)


width_pw = 0.13
height_pw = 2.15

# Create the rectangle patch
#rectangle_pw = mpatches.Rectangle((33.89, -0.17), width_pw, height_pw,
#                               facecolor='thistle', edgecolor='purple', linewidth=2, alpha=0.7)

# Add the rectangle to the axes
#axes[1].add_patch(rectangle_pw)

t400 = False
t300 = False
t500 = False
t150 = False
t125 = False

for i in range(8):
    dist = distances[i]
    #color = colormap_oxy(norm(dist)) # Assigns each cast a color in the spectrum
    color = colormap_turb(norm_turb(turb[i]))
    #color = colormap_oxy(norm_depth(-depth[i]))
    #color = "darkgray"
    if i in [10, 11, 12, 18, 22, 23, 24]:
        color = "red"
    axes[1].scatter(SA[i], CT[i], c=color, s=3, zorder=5)
    
    for j in range(len(depth[i])):
        if float(depth[i][j])==500: 
            #axes[1].scatter(SA[i][j], CT[i][j], color ="deeppink", s = 10, zorder = 3)
            if t500 == False: 
                #axes[1].text(SA[i][j]+0.05, CT[i][j]-0.2, "500m", color = "deeppink")
                t500=True
    
    # for j in range(len(depth[i])):
    #     if float(depth[i][j])==400: 
    #         axes[1].scatter(SA[i][j], CT[i][j], color ="magenta", s = 10, zorder = 3)
    #         if t400 == False: 
    #             axes[1].text(SA[i][j]+0.03, CT[i][j]-0.2, "400m", color = "magenta")
    #             t400=True
                
    for j in range(len(depth[i])):
        if float(depth[i][j])==300: 
            #axes[1].scatter(SA[i][j], CT[i][j], color ="blueviolet", s = 10, zorder = 3)
            if t300 == False: 
                #axes[1].text(SA[i][j]+0.05, CT[i][j]-0.4, "300m", color = "blueviolet")
                t300=True
                
    
    for j in range(len(depth[i])):
        if float(depth[i][j])==150: 
            #axes[1].scatter(SA[i][j], CT[i][j], color ="blue", s = 10, zorder = 3)
            if t150 == False: 
                #axes[1].text(SA[i][j]-0.2, CT[i][j]+0.4, "150m", color = "blue")
                t150=True
                
    # for j in range(len(depth[i])):
    #     if float(depth[i][j])==125: 
    #         axes[1].scatter(SA[i][j], CT[i][j], color ="cyan", s = 10, zorder = 3)
    #         if t125 == False: 
    #             axes[1].text(SA[i][j]-0.2, CT[i][j]+0.4, "125m", color = "cyan")
    #             t125=True
    
    
axes[1].plot(SA[8], CT[8], c='k', label="Cast 9 (shelf)", zorder=5)

axes[1].set_xlim(33.5, 34.9)
axes[1].set_ylim(-0.5, 4)
axes[1].scatter(34.237632875, 1.8995849875, marker = "*", s=70, color = "red", zorder = 6, label = "Turbidity Peaks")
axes[1].scatter(34.4352518, 2.0476521400000003, marker = "*", s=70, zorder = 6, color = "red")

axes[1].plot(sx, runoff_line(sx), color ="orange", linestyle = "dashed")
axes[1].plot(sx, melting_line(sx), color ="pink", linestyle = "dashed")

axes[1].legend(loc="lower right")
#axes[1].set_title(r"$\mathbf{Fjord:}$" + "\nColored by Turbidity")
#axes[1].set_title(r"$\mathbf{Fjord:}$")
axes[1].set_ylabel("Conservative Temperature [°C]")
axes[1].set_xlabel("Absolute Salinity [g/kg]")


cbar_turb = fig.colorbar(sm_turb, ax=axes[1], orientation='vertical')
#cbar_turb = fig.colorbar(sm_depth, ax=axes[1], orientation='vertical')
#cbar_turb = fig.colorbar(sm, ax=axes[1], orientation='vertical')
#cbar_turb.set_label("Turbidity [NTU]")
cbar_turb.set_label("Turbidity [NTU]")



plt.tight_layout()
plt.show()

