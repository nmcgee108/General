#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 14:09:26 2025

@author: nataliemcgee
"""

import xarray as xr
import numpy as np
from netCDF4 import Dataset
from pyproj import Geod
import pandas as pd
import gsw
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from matplotlib.colors import PowerNorm

plt.rcParams['font.size'] = 20


# Load the CSV file into a pandas DataFrame
bath_df = pd.read_csv('/Users/nataliemcgee/Documents/Upernavik Data/Bathymetry Data/new_fjord_bathymetry1.csv')


# Extract columns for plotting
section_distances = bath_df['Distance_km']  
bed = bath_df['Bed_Elevation_m'] 
section_lats = bath_df['Latitude']  
section_lons = bath_df['Longitude'] 

# Add CTD lines
geod = Geod(ellps='WGS84')

start_cast = 1
end_cast = 9

# Ensure cast numbers and lat/lon arrays align correctly
cast_numbers = np.arange(start_cast, end_cast + 1)  # [11, 12, ..., 19]
n_casts = len(cast_numbers)

#CTD locations

#fjord
ctd_lats = [
    72.91,
    72.8849,
    72.88805,
    72.8822,
    72.91523,
    72.9404,
    73.00466,
    73.017833,
    73.127]


ctd_lons = [
    -55.01,
    -55.0707,
    -55.0057,
    -54.9242,
    -55.454,
    -55.6726,
    -56.08366,
    -56.330163,
    -57.3294]

# #Trough
# ctd_lats = [
#     72.7735,
#     72.83995,
#     72.91643,
#     72.9867,
#     73.05838,
#     73.12492,
#     73.19801,
#     73.27043,
#     73.3421]

# ctd_lons = [
#     -58.58975,
#     -58.7702,
#     -58.92601,
#     -59.146,
#     -59.34056,
#     -59.504433,
#     -59.70404,
#     -59.90671,
#     -60.076983]

def find_closest_points(ctd_lats, ctd_lons, section_lats, section_lons, section_distances):
    ctd_dist_along_section = []

    for i in range(len(ctd_lats)):
        # Compute distances to all section points from the i-th CTD
        _, _, dists = geod.inv(
            np.full_like(section_lons, ctd_lons[i]), #make an array shaped like section_lons but fill with one value of CTD-lons
            np.full_like(section_lats, ctd_lats[i]),
            section_lons, #compare item by item with section_lons itself
            section_lats)
        
        closest_point = np.argmin(np.abs(dists))
        ctd_dist_along_section.append(section_distances[closest_point]) 


    return np.array(ctd_dist_along_section)

ctd_dist_along_section=find_closest_points(ctd_lats, ctd_lons, section_lats, section_lons, section_distances) #all distances starting with cast 1 to cast 8 in cast order


uc_patch_file = "/Users/nataliemcgee/Documents/Upernavik Data/Final CTD Datasets/uc_patch_dataset_new.nc"


with Dataset(uc_patch_file, 'r') as ncfile:

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
    

# depth_bin : list of 1D arrays of depths, one per cast
# SA, CT, theta, turb, oxy, etc : list of 1D arrays per cast 
# latitudes, longitudes, ctd_distances : arrays of shape (n_casts,)

# Build the full common depth list
all_depths = np.unique(np.concatenate(depth_bin))
all_depths = np.sort(all_depths)
n_depths = len(all_depths)

# Create NaN-filled arrays
SA_aligned     = np.full((n_casts, n_depths), np.nan)
CT_aligned     = np.full((n_casts, n_depths), np.nan)
theta_aligned  = np.full((n_casts, n_depths), np.nan)
turb_aligned   = np.full((n_casts, n_depths), np.nan)
oxy_aligned    = np.full((n_casts, n_depths), np.nan)
fluor_aligned  = np.full((n_casts, n_depths), np.nan)
sigma0_aligned = np.full((n_casts, n_depths), np.nan)

max_depths = {} # The bottom depth of each cast

for cast in range(start_cast, end_cast + 1):
    i = cast - 1  # zero-based index
    cast_depths = np.array(depth_bin[i])
    max_depths[cast] = np.nanmax(cast_depths)
    
print("Max depths by cast:")
for k in sorted(max_depths):
    print(k, max_depths[k])
    
# Populate aligned arrays by matching depth values
for idx, i in enumerate(range(start_cast-1, end_cast)):
    if i == 9 :
        print("cast 10 skipped")
        continue
    cast_depths = np.array(depth_bin[i])  # All depth bins in this cast
    #max_depths.append(max(cast_depths))
    #print(max_depths)
    indices = np.searchsorted(all_depths, cast_depths) # Finds indices where depth bins align with all_depth array

    SA_i     = np.asarray(SA[i])
    CT_i     = np.asarray(CT[i])
    theta_i  = np.asarray(theta[i])
    turb_i   = np.asarray(turb[i])
    oxy_i    = np.asarray(oxy[i])
    fluor_i  = np.asarray(fluor[i])

    #SA_aligned[idx, indices]     = savgol_filter(SA_i, window_length=15, polyorder=2)
    SA_aligned[idx, indices]     = SA_i
    CT_aligned[idx, indices]     = savgol_filter(CT_i, window_length=15, polyorder=2)
    theta_aligned[idx, indices]  = savgol_filter(theta_i, window_length=15, polyorder=2)
    turb_aligned[idx, indices]   = savgol_filter(turb_i, window_length=10, polyorder=2)
    oxy_aligned[idx, indices]    = savgol_filter(oxy_i, window_length=15, polyorder=2)
    fluor_aligned[idx, indices]  = savgol_filter(fluor_i, window_length=15, polyorder=2)
    
    # Compute sigma0 and smooth vertically
    valid = ~np.isnan(SA_i) & ~np.isnan(CT_i)
    sigma0_vals = gsw.sigma0(SA_i[valid], CT_i[valid])
    smoothed_sigma0 = savgol_filter(sigma0_vals, window_length=100, polyorder=3)
    sigma0_aligned[idx, indices[valid]] = smoothed_sigma0

        
# Build the xarray Dataset
section = xr.Dataset(
    data_vars={
        "Absolute Salinity": (["cast", "depth"], SA_aligned),
        "Conservative Temperature": (["cast", "depth"], CT_aligned),
        "Potential Temperature": (["cast", "depth"], theta_aligned),
        "Turbidity": (["cast", "depth"], turb_aligned),
        "Oxygen": (["cast", "depth"], oxy_aligned),
        "Fluorescence": (["cast", "depth"], fluor_aligned),
        "Potential Density Anomaly": (["cast", "depth"], sigma0_aligned)},
    coords={
        "cast": cast_numbers,
        "depth": all_depths,
        "lat": ("cast", ctd_lats),
        "lon": ("cast", ctd_lons),
        "distance": ("cast", ctd_dist_along_section)},
    attrs={
        "description": "Aligned CTD cast data from Upernavik Fjord",
        "note": "Variables aligned to shared depth grid; NaNs indicate missing data"})



def fill_gaps(ds, max_gap=4):
    """
    Fill vertical gaps (up to `max_gap` meters) in each cast using linear interpolation.
    
    Parameters:
    - ds: xarray.Dataset with dims (cast, depth)
    - max_gap: maximum gap size (in meters) to fill
    
    Returns:
    - xarray.Dataset with small gaps filled
    """
    filled_ds = ds.copy()

    for var in ds.data_vars:
        if ds[var].dims != ("cast", "depth"):
            continue  # Skip non-profile variables

        data = ds[var].copy()

        for i in range(data.sizes["cast"]):
            profile = data.isel(cast=i).values
            valid = ~np.isnan(profile)

            if valid.sum() < 2:
                continue  # Not enough data to interpolate

            # Find runs of NaNs
            is_nan = ~valid
            gap_start = None
            for j in range(1, len(profile) - 1):
                if is_nan[j] and not is_nan[j - 1]:
                    gap_start = j
                if gap_start is not None and not is_nan[j] and is_nan[j - 1]:
                    gap_end = j
                    gap_len = gap_end - gap_start
                    if gap_len <= max_gap:
                        # Interpolate over the gap
                        y0, y1 = profile[gap_start - 1], profile[gap_end]
                        interp = np.linspace(y0, y1, gap_len + 2)[1:-1]
                        profile[gap_start:gap_end] = interp
                    gap_start = None  # Reset after closing a gap

            data[i, :] = profile

        filled_ds[var] = data

    return filled_ds

section = fill_gaps(section)


def plot_section(section_dataset, var_name, label, 
                 casts_to_plot=None, reference_cast=None, vmin=None, vmax=None,
                 cmap='viridis', 
                 levels=100, 
                 gamma =0.5,
                 sigma_levels=[26, 26.5, 27, 27.5]):
    """
    Plot a 2D section of a variable with isopycnals and cast markers.
    
    Parameters:
    - section: xarray.Dataset with dims (cast, depth)
    - var_name: variable to plot (e.g., 'turb', 'CT', 'SA')
    - label: label for the colorbar
    - casts_to_plot: list of cast indices to include (default = all)
    - cmap: matplotlib colormap
    - levels: contour levels (int or list)
    - sigma_levels: list of sigma0 contour levels
    """
    # Subset dataset if specific casts requested
    if casts_to_plot is not None:
        section = section_dataset.sel(cast=casts_to_plot)
        section = section.sortby("distance") #put them in the correct order

    # Extract variables
    var = section[var_name]
    sigma = section["Potential Density Anomaly"]
    dist = section["distance"]
    depth = section["depth"]
    
    if vmin == None:
        vmin = float(var.min())

    if vmax == None: 
        vmax = float(var.max())

    # Build meshgrid
    X, Y = np.meshgrid(dist, -depth)

    fig, ax = plt.subplots(figsize=(12, 6))
    cf = ax.contourf(X, Y, var.T, levels=levels, norm=PowerNorm(gamma=gamma), cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(cf, ax=ax, format="%.2f")
    cbar.set_label(label, labelpad=15)

    # Overlay sigma0 contours
    cs = ax.contour(X, Y, sigma.T, levels=sigma_levels, colors='white', linewidths=1)
    ax.clabel(cs, inline=True, fmt="%.2f", inline_spacing = 70)

    # Add cast vertical lines at actual max sample depth
    for i, cast in enumerate(section.cast.values):
        ax.vlines(dist.values[i], 0, -max_depths[int(cast)], color='k', linestyle='-')
        plt.scatter(dist.values[i], 15, marker = "^", color = "red")
        plt.text(dist.values[i]+1, 40, cast)
    
    plt.tick_params(axis='both', which='major')
    
    
    
    
    

    #########
    # Add shelf cast for reference
    if reference_cast != None:
        cast9 = section_dataset.isel(cast=8)
        cast9 = fill_gaps(cast9)
    
        # Extract variable and depth
        var_profile = cast9[var_name].values  # shape: (depth,)
        sigma_profile = cast9["Potential Density Anomaly"].values
        depth_ref = cast9["depth"].values
        
        # Build 2-column "block" for plotting (same data in both columns)
        var_ref = np.column_stack([var_profile, var_profile])  # shape: (depth, 2)
        sigma_ref = np.column_stack([sigma_profile, sigma_profile])
        
        # Choose position offset from main section
        dist_center = float(cast9["distance"].values)
        dist_ref = np.array([dist_center - 5, dist_center + 5])  # two fake "x" positions
        
        # Build meshgrid
        X_ref, Y_ref = np.meshgrid(dist_ref, -depth_ref)
        
        # Plot filled contour
        #ax.contourf(X_ref, Y_ref, var_ref, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.contourf(X_ref, Y_ref, var_ref, levels=levels, norm=PowerNorm(gamma=gamma), cmap=cmap, vmin=vmin, vmax=vmax)
        
        # overlay density contours
        ax.contour(X_ref, Y_ref, sigma_ref, levels=sigma_levels, colors='white', linewidths=1)
        #ax.clabel(cs_ref, inline=True, fmt="%.1f", fontsize=9, inline_spacing=10)
    
        ax.vlines(dist_center, 0, -max_depths[9], color='k', linestyle='-')
        plt.scatter(dist_center, 15, marker = "^", color = "red")
        plt.text(dist_center+1, 40, "9")


    ax.plot(section_distances, bed, color="gray")
    ax.fill_between(section_distances, bed, min(bed) - 50, color="gray", zorder=0)
    ax.set_xlim(-0.2, 110)
    ax.invert_xaxis()

    plt.text(77, -1020, "Mouth")
    plt.text(15, -1020, "Glacier")
    plt.text(108, -700, "Shelf")
    ax.set_ylim(-1050, 150)
    ax.set_xlabel("Distance from Glacier [km]")
    ax.set_ylabel("Depth [m]")
    #ax.set_title(f"{var_name} Section Upernavik Fjord", fontsize=20)
    plt.tight_layout()
    plt.show()



plot_section(section, "Absolute Salinity", "Absolute Salinity [g/kg]", 
                 casts_to_plot=np.arange(2,9), # Cast num
                 reference_cast=9,
                 cmap='viridis', 
                 levels=200, 
                 vmin = 29,
                 vmax = 34.5,
                 gamma = 1,
                 sigma_levels=[26.5, 27, 27.5])

plot_section(section, "Turbidity", "Turbidity [NTU]", 
                 casts_to_plot=np.arange(2,9), # Cast num
                 reference_cast=9,
                 cmap='viridis', 
                 levels=100, 
                 gamma = 0.7,
                 sigma_levels=[26.5, 27, 27.5])

plot_section(section, "Oxygen", "Oxygen [ml/l]", 
                 casts_to_plot=np.arange(2,9), # Cast num
                 reference_cast=9,
                 cmap='viridis', 
                 levels=200, 
                 vmin = 5,
                 vmax = 12.2,
                 gamma = 0.45,
                 sigma_levels=[26.5, 27, 27.5])

plot_section(section, "Conservative Temperature", "Conservative Temperature [°C]", 
                 casts_to_plot=np.arange(2,9), # Cast num
                 reference_cast=9,
                 cmap='viridis', 
                 levels=100, 
                 vmin = 0,
                 vmax = 3.7,
                 gamma = 0.7,
                 sigma_levels=[26.5, 27, 27.5])#[26.5, 27, 27.1, 27.2, 27.3, 27.4, 27.5, 27.55])

plot_section(section, "Fluorescence", r"Fluorescence [mg/m$^3$]", 
                 casts_to_plot=np.arange(2,9), # Cast num
                 reference_cast=9,
                 cmap='viridis', 
                 levels=100, 
                 gamma = 0.5,
                 sigma_levels=[26.5, 27, 27.5])
