#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot a bathymetric section made of two line segments with data from bedmachinegreenland 
and create a csv file for future use

Created on Mon Jun 16 17:16:03 2025

@author: nataliemcgee
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pyproj
from pyproj import Geod



ds = xr.open_dataset("/Users/nataliemcgee/Documents/Upernavik Project/Upernavik Data/BedMachineGreenland-v5.nc")
bed = ds['bed'].values
x = ds['x'].values
y = ds['y'].values

# Create 2D meshgrid
xx, yy = np.meshgrid(x, y)


# Define projections
stere = pyproj.CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +datum=WGS84")
geodetic = pyproj.CRS("EPSG:4326")  # WGS84 lat/lon

# Create transformers
to_geodetic = pyproj.Transformer.from_crs(stere, geodetic, always_xy=True)
to_projected = pyproj.Transformer.from_crs(geodetic, stere, always_xy=True)


# Define start and end coordinates
start_lat, start_lon = 72.89, -54.7
mid1_lat, mid1_lon = 72.89061111111111, -55.23407777777778
end_lat, end_lon = 73.24633333333334, -58.31155555555556

# Convert start/end to projected x/y
x1, y1 = to_projected.transform(start_lon, start_lat)
x2, y2 = to_projected.transform(mid1_lon, mid1_lat)
#x3, y3 = to_projected.transform(mid2_lon, mid2_lat)
x3, y3 = to_projected.transform(end_lon, end_lat)


# Sample points along the line
N = 500
x_line1 = np.linspace(x1, x2, N)
y_line1 = np.linspace(y1, y2, N)

x_line2 = np.linspace(x2, x3, N)
y_line2 = np.linspace(y2, y3, N)

#x_line3 = np.linspace(x3, x4, N)
#y_line3 = np.linspace(y3, y4, N)

# Interpolate bed elevation along line
bed_interp = RegularGridInterpolator((y, x), bed, bounds_error=False, fill_value=np.nan)

bed_along_line1 = bed_interp(np.array([y_line1, x_line1]).T)
bed_along_line2 = bed_interp(np.array([y_line2, x_line2]).T)
#bed_along_line3 = bed_interp(np.array([y_line3, x_line3]).T)

#bed_along_full_line = bed_along_line1
bed_along_full_line = np.concatenate((bed_along_line1, bed_along_line2))#, bed_along_line3)) # Combine bed data

# Compute distance along the line
lons1, lats1 = to_geodetic.transform(x_line1, y_line1)
lons2, lats2 = to_geodetic.transform(x_line2, y_line2)
#lons3, lats3 = to_geodetic.transform(x_line3, y_line3)
geod = Geod(ellps='WGS84')

#lons,lats = lons1, lats1
lons, lats = np.concatenate((lons1, lons2)), np.concatenate((lats1, lats2)) # Combine lats and lons

distances = [0]
#for i in range(1, 3*N):
for i in range(1, 2*N):
    _, _, d = geod.inv(lons[i-1], lats[i-1], lons[i], lats[i])
    distances.append(distances[-1] + d / 1000)  # Convert to km
    
bed_array = np.column_stack((distances, bed_along_full_line, lats, lons))
np.savetxt('new_fjord_bathymetry1.csv', bed_array, delimiter=',', header='Distance_km,Bed_Elevation_m,Latitude,Longitude',
    comments='')
print("file saved. format is distance, depth, lat, lon")

    
#---------------------------------------------------------
# now add the ctd stations
#---------------------------------------------------------

ctd_lats = [72.91, 72.8849, 72.88805, 72.8822, 72.91523, 72.9404, 73.00466, 73.017833]
ctd_lons = [-55.01, -55.0707, -55.0057, -54.9242, -55.454, -55.6726, -56.08366, -56.330163]

# #ctd_lats = [
#     72.7735,
#     72.83995,
#     72.91643,
#     72.9867,
#     73.05838,
#     73.12492,
#     73.19801,
#     73.27043,
#     73.3421]

# #ctd_lons = [
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
            np.full_like(section_lons, ctd_lons[i]),
            np.full_like(section_lats, ctd_lats[i]),
            section_lons,
            section_lats)
        
        closest_point = np.argmin(dists)
        ctd_dist_along_section.append(section_distances[closest_point])

    return np.array(ctd_dist_along_section)

ctd_dist_along_section=find_closest_points(ctd_lats, ctd_lons, lats, lons, distances) #all distances starting with cast 1 to cast 8 in cast order



# Plot the cross-section
plt.figure(figsize=(10, 5))
plt.plot(distances, bed_along_full_line, color='steelblue')
plt.vlines(ctd_dist_along_section, -800, 0, color = 'k')
plt.scatter(ctd_dist_along_section, np.zeros_like(ctd_dist_along_section), marker = "^", color = "red")
plt.fill_between(distances, bed_along_full_line, min(bed_along_full_line)-50) #shade below the bed
plt.gca().invert_xaxis()
plt.xlabel("Distance (km)") 
plt.ylabel("Bed elevation (m)")
plt.title("Bathymetric Section")
plt.tight_layout()
plt.show()
