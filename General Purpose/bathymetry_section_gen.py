#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GENERAL PUROPSE

Plot a bathymetric section made of a line segment (lat, lon) with data from BedMachineGreenland 
and save bed depths along the line as a .csv file for future use.

BedMachine files expected in NetCDF format and North polar stereographic projection (default).
This projection contains data in an x, y grid that is not lat, lon and must be transformed.

Recommended environment: upernavik_env.yml

Created on Mon Jun 16 17:16:03 2025

@author: nataliemcgee
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pyproj
from pyproj import Geod

# Load the BedMachine NetCDF dataset
ds = xr.open_dataset("/Users/nataliemcgee/Downloads/Upernavik Project/Upernavik Data/BedMachineGreenland-v5.nc") 
bed = ds['bed'].values # Extract bed depth values
x = ds['x'].values  # Extract x location (from projection)
y = ds['y'].values # Extract y location (from projection)

# Create 2D meshgrid
xx, yy = np.meshgrid(x, y) 



# Define the North Polar Stereographic projection of the BedMachine data. 
# Change this if your data uses another projection
stere = pyproj.CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +datum=WGS84") 

# Define the lat/lon (WGS84 geodetic) projection
geodetic = pyproj.CRS("EPSG:4326")  




# Create transformers
to_geodetic = pyproj.Transformer.from_crs(stere, geodetic, always_xy=True)
to_projected = pyproj.Transformer.from_crs(geodetic, stere, always_xy=True)


# Define start and end coordinates of section line
start_lat, start_lon = 73.055, -56.546 # Negative values indicate West lons
end_lat, end_lon = 72.89, -54.48

# Convert start/end to projected x/y to map onto BedMachine
x1, y1 = to_projected.transform(start_lon, start_lat)
x2, y2 = to_projected.transform(end_lon, end_lat)


# Sample points along the line
N = 500
x_line = np.linspace(x1, x2, N) # Array of x-values of the line in stereo projection
y_line = np.linspace(y1, y2, N)


# Interpolate bed elevation along line
bed_interp = RegularGridInterpolator((y, x), bed, bounds_error=False, fill_value=np.nan)
bed_along_line = bed_interp(np.array([y_line, x_line]).T)

# Compute distance along the line for plotting purposes
lons, lats = to_geodetic.transform(x_line, y_line) 

# Define the geod and desired ellipsoid. WGS84 is standard for Earth
# Used to calculate great circle distances between lat, lon points 
geod = Geod(ellps='WGS84')

distances = [0]

# Calculate distance between every point and sum to find total distance from starting point
for i in range(1, N):
    _, _, d = geod.inv(lons[i-1], lats[i-1], lons[i], lats[i]) # geod.inv returns three arguments; only the third (distance, m) is useful 
    distances.append(distances[-1] + d / 1000)  # Convert to km
    
    
# Plot the cross-section
plt.figure(figsize=(10, 5))
plt.plot(distances, bed_along_line, color='steelblue')
plt.fill_between(distances, bed_along_line, min(bed_along_line)-50) # Shade below the bed
plt.xlabel("Distance along section (km)") 
plt.ylabel("Bed elevation (m)")
plt.title("Bathymetric Section")
plt.tight_layout()
plt.show()
    
    
#---------------------------------------------------------
# Save depth, distance as a .csv file
#---------------------------------------------------------  

filename = 'fjord_section_bathymetry.csv'
bed_array = np.column_stack((distances, bed_along_line))
np.savetxt(filename, bed_array, delimiter=',', header='Distance_km,Bed_Elevation_m',
    comments='')

print(f"File saved at {filename}. Format is distance, depth")


