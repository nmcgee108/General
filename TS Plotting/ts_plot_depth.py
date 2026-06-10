#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 15:57:57 2025

@author: nataliemcgee
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import numpy as np
import gsw

uc_patch_file = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Final CTD Datasets/uc_patch_dataset_new.nc"

ncfile = Dataset(uc_patch_file, 'r')

#print all the variables and their dimensions:
#print(ncfile)

#print the names (keys) of the variables
#print(ncfile.variables.keys())

#upernavik fjord are casts 1-8, possibly 9 which is on the shelf outside the mouth

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
norm_depth = mcolors.Normalize(vmin=-750, vmax=0)
norm_castnum = mcolors.Normalize(vmin=10, vmax=18)

sm_oxy = plt.cm.ScalarMappable(cmap=colormap_oxy, norm=norm_oxy)
sm_turb = plt.cm.ScalarMappable(cmap=colormap_turb, norm=norm_turb)
sm_depth = plt.cm.ScalarMappable(cmap=colormap_oxy, norm=norm_depth)
sm_castnum = plt.cm.ScalarMappable(cmap=colormap_oxy, norm=norm_castnum)

minS = 28.5
maxS = 35
minT = -2
maxT = 7.4

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

# Create figure and custom layout: 2 plots + 1 for colorbar
fig = plt.figure(figsize=(12, 5))  # Slightly wider to accommodate colorbar
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 0.05], wspace=0.3)

# Create subplots
axes = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]

########### shelf plot ########
contour = axes[0].contour(S, T, PDEN, levels=np.arange(23, 28), colors='grey')
axes[0].clabel(contour, inline=True, fmt='%1.0f') 

for i in range(11, 19):
    color = colormap_oxy(norm_depth(-depth[i]))
    axes[0].scatter(SA[i], CT[i], c=color, s=2, marker = ".")
    #axes[0].plot(SA[i], CT[i], c=color)
 
color = colormap_oxy(norm_depth(-depth[30]))
axes[0].scatter(SA[30], CT[30], c=color, s=2, marker = "^", label="Cast 31 (shelf)")
color = colormap_oxy(norm_depth(-depth[31]))
axes[0].scatter(SA[31], CT[31], c=color, s=10, marker = "^", label="Cast 32 (shelf)")

#axes[0].plot(SA[12], CT[12], c='k', label="Cast 9 (shelf)")
axes[0].plot(sx, runoff_line(sx), color="orange", linestyle="dashed", label="Runoff Line")
axes[0].plot(sx, melting_line(sx), color="pink", linestyle="dashed", label="Melting Line")

axes[0].set_title("Shelf and Trough")
axes[0].set_xlabel("Absolute Salinity [g/kg]")
axes[0].set_ylabel("Conservative Temperature [°C]")
axes[0].set_xlim(33.5, 35)
axes[0].set_ylim(0, 5)
axes[0].legend()

######### fjord plot ########
contour = axes[1].contour(S, T, PDEN, levels=np.arange(23, 28), colors='grey')
axes[1].clabel(contour, inline=True, fmt='%1.0f')

for i in range(8):
    color = colormap_oxy(norm_depth(-depth[i]))
    axes[1].scatter(SA[i], CT[i], c=color, s=2, marker=".")

axes[1].plot(SA[8], CT[8], c='k')
axes[1].plot(SA[51], CT[51], c='k', label="Cast 9 (shelf)")
axes[1].scatter(34.237632875, 1.8995849875, marker="*", s=20, color="red", label="Turbidity Peaks")
axes[1].scatter(34.4352518, 2.0476521400000003, marker="*", s=20, color="red")

axes[1].plot(sx, runoff_line(sx), color="orange", linestyle="dashed")
axes[1].plot(sx, melting_line(sx), color="pink", linestyle="dashed")

axes[1].set_title("Fjord")
axes[1].set_xlabel("Absolute Salinity [g/kg]")
axes[1].set_ylabel("Conservative Temperature [°C]")
axes[1].set_xlim(32.5, 35)
axes[1].set_ylim(-2, 5)
axes[1].legend(loc="lower right")

########### Shared Colorbar ########
# Create a dedicated axis for the colorbar
cax = fig.add_subplot(gs[2])
cbar_turb = fig.colorbar(sm_depth, cax=cax, orientation='vertical')
cbar_turb.set_label("Depth [m]")

fig.suptitle("Upernavik Fjord", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.95])  # Leave space for suptitle
plt.show()

