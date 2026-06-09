#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 15:49:24 2025

@author: nataliemcgee
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import gsw
from pyproj import Geod
from full_plume import run_plume
from analytical_plume_gen import analytical_plume

# Load the CTD NetCDF dataset
ctd_file = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/Final CTD Datasets/uc_patch_dataset_new.nc"

with Dataset(ctd_file, 'r') as ncfile:
    
    # Extract cast and depth-indexed variables (many values per cast)
    # Creates a 2D numpy array of values indexed by cast number and depth bin (missing data filled with NaN)
    lat = np.ma.filled(ncfile.variables['Latitude'][:], np.nan) 
    lon = np.ma.filled(ncfile.variables['Longitude'][:], np.nan)
    depth = np.ma.filled(ncfile.variables['Depth Bin'][:], np.nan)
    SA = np.ma.filled(ncfile.variables['Absolute Salinity'][:], np.nan)
    CT = np.ma.filled(ncfile.variables['Conservative Temp'][:], np.nan)
    
#-----------------------------
# THINGS TO ENTER
#-----------------------------
cast_num = 3
plume_depth = 659
Q_discharge = 60 # subglacial discharge (m3/s)
width = 500 # (m)

#-----------------------------
# Full plume model
#-----------------------------
Q0 = Q_discharge/width # (m2/s)
depth_mask = np.where(depth[cast_num-1]<plume_depth)
zi = -depth[cast_num-1][depth_mask][::-1] # Must input values from deepest to shallowest!
xi = np.zeros_like(zi)
Ta = CT[cast_num-1][depth_mask][::-1]   # ambient fjord temperature at zi
Sa = SA[cast_num-1][depth_mask][::-1] # ambient fjord salinity at zi
# Sa=np.average(Sa[:5])*np.ones_like(zi) # Ambient properties: average the first ?m above the plume depth
# Ta=np.average(Ta[:5])*np.ones_like(zi)
Na = np.ones_like(zi)*10 # ambient fjord nitrate at zi
alpha = 0.1 # entrainment coefficient

if max(depth[cast_num-1])<plume_depth:
    print("****Warning: Chosen plume depth exceeds CTD data! T, S extrapolated at depth*****")
    zi = np.insert(zi, 0, -plume_depth)
    xi = np.insert(xi, 0, 0)
    Ta = np.insert(Ta, 0, np.average(Ta[:10]))
    Sa = np.insert(Sa, 0, np.average(Sa[:10]))
    Na = np.insert(Na, 0, 10)

sol=run_plume(zi, xi, Ta, Sa, Na, Q0, alpha)

#-----------------------------
# Analytical plume model
#-----------------------------

S_AW=np.average(Sa[:50]) # Ambient properties: average the first ?m above the plume depth
T_AW=np.average(Ta[:50])
Tplume, Splume, AWp, SGDp, SMWp, Q_AW, Q_SMW = analytical_plume(T_AW=T_AW,
                                                                S_AW=S_AW,
                                                                Q_SGD=Q_discharge,
                                                                h_gl=plume_depth,
                                                                w=width)

#---------------------------------------------------------
# Set up color coding by distance
#---------------------------------------------------------

# Define the geod and desired ellipsoid. WGS84 is standard for Earth
# Used to calculate great circle distances between lat, lon points 
geod = Geod(ellps='WGS84')

lon_ref, lat_ref = 54.4036 , 72.9435
distances = [] # Distances will fill the list indexed by cast order

# Calculate distance between every point and reference point
for i in range(8):
    _, _, d = geod.inv(lon_ref, lat_ref, lon[i], lat[i]) # geod.inv returns three arguments; only the third (distance, m) is useful 
    distances.append(d / 1000)  # Convert to km


max_dist = max(distances) 
norm_distance = mcolors.Normalize(vmin=0, vmax=max_dist)  # Normalize the colormap with a max and min value

#---------------------------------------------------------
#---------------------------------------------------------

norm_depth = mcolors.Normalize(vmin=-750, vmax=0)

colormap = plt.colormaps['viridis']
sm_depth = plt.cm.ScalarMappable(cmap=colormap, norm=norm_depth)
sm_distance = plt.cm.ScalarMappable(cmap=colormap, norm=norm_distance)

minS = 32.5
maxS = 35
minT = -2
maxT = 4

sx = np.arange(minS, maxS + 0.1, 0.1)
ty = np.arange(minT, maxT + 0.1, 0.1)


# Create grid of salinity and temperature
S, T = np.meshgrid(sx, ty)
PDEN = gsw.rho(S, T, 0) - 1000  # Potential density anomaly at reference pressure = 0

#mixing line params
T_deepwater = T_AW #2.550228
S_deepwater = S_AW #34.697933
T_ice_eff = -90
    
def runoff_line(x):
    slope = T_deepwater/S_deepwater
    return slope*x

def melting_line(x):
    slope = (T_deepwater-T_ice_eff)/S_deepwater
    return slope*x+T_ice_eff

# Create figure and custom layout: 2 plots + 1 for colorbar
fig, axes = plt.subplots(1, 1, figsize=(7, 6))  # Slightly wider to accommodate colorbar


######### fjord plot ########

contour = axes.contour(S, T, PDEN, levels=np.arange(26, 28, 0.5), colors='grey')
axes.clabel(contour, inline=True, fmt='%0.1f')

axes.plot(sx, runoff_line(sx), color="orange", linestyle="dashed", zorder = 0)
axes.plot(sx, melting_line(sx), color="pink", linestyle="dashed", zorder = 0)

for i in range(8):
    if i == 8: axes.plot(SA[i], CT[i], color ='k')
    else:
        color = colormap(norm_distance(distances[i]))
        axes.scatter(SA[i], CT[i], color=color, s=2, marker=".")
    
axes.scatter(S_AW, T_AW, marker = "o", s = 30, color = "k", label = "Calculated Ambient")

#axes.plot(sol["S"], sol["T"], c = "darkblue", linestyle = "dotted")
axes.scatter(sol["S"][-1], sol["T"][-1], c="darkblue", marker = "s", s=60, label = "Full Model Output")

#axes.plot(SA[8], CT[8], c='k')
#axes.plot(SA[10], CT[10], c='k')
#axes.plot(SA[11], CT[11], c='k')
#axes.plot(SA[12], CT[12], c='k')
#axes.plot(SA[13], CT[13], c='k')
#axes.plot(SA[14], CT[14], c='k')

S_FPW = 34.13 #33.6 #33.5   # Fjord Polar Water (not shelf PW. See Muilwijk definitions)
T_FPW = 1.6 #0.5 #1.8
P_FPW = 0.41



filled = axes.contourf(
    S, T, PDEN,
    levels=[26.98, 27.43],  # Just the range you want shaded
    colors=['lightskyblue'],  # Or use a colormap like cmap='Greys'
    alpha=0.3,             # Optional: make it slightly transparent
    extend='neither',       # Controls if arrows show on the ends
    zorder =0
    )

# filled = axes.contourf(
#     S, T, PDEN,
#     levels=[26.95, 27.2],
#     colors=['pink'],  # Or use a colormap like cmap='Greys'
#     alpha=0.3,             # Optional: make it slightly transparent
#     extend='neither' ,      # Controls if arrows show on the ends
#     zorder=0
#     )

axes.scatter(Splume, Tplume, marker = "s", s = 60, color = "darkviolet", label = "Analytic Model Output", zorder = 4)
axes.scatter(Splume*(1-P_FPW)+ S_FPW*P_FPW, Tplume*(1-P_FPW) + T_FPW*P_FPW, marker = "^", s = 80, color = "darkviolet", label = "Mod. Analytic Output", zorder =4)

axes.scatter(34.237632875, 1.8995849875, marker="*", s=40, color="red")
axes.scatter(34.4352518, 2.0476521400000003, marker="*", s=40, color="red")

#axes.scatter(34.13, 1.6, marker="*", s=40, color="red")

axes.set_title(f"Plume Depth = {plume_depth} m")
axes.set_xlabel("Absolute Salinity [g/kg]")
axes.set_ylabel("Conservative Temperature [°C]")
axes.set_xlim(33.5, 35)
axes.set_ylim(-0.2, 3.5)
axes.legend(loc="lower right")

########### Shared Colorbar ########
# Create a dedicated axis for the colorbar
cbar_turb = fig.colorbar(sm_distance, ax = axes, orientation='vertical')
cbar_turb.set_label("Distance from Glacier [km]")

plt.show()

