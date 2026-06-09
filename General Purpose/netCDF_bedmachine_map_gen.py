#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GENERAL PURPOSE

Create a map with data from BedMachineGreenland and plot points in lat/lon.

BedMachine files expected in NetCDF format and North polar stereographic projection (default).
This projection contains data in an x, y grid that is not lat, lon and must be transformed.

Recommended environment: upernavik_env.yml

Created on Mon Jun 16 12:27:38 2025

@author: nataliemcgee
"""

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pyproj
import numpy as np

# Load the BedMachine NetCDF dataset
ds = xr.open_dataset("/Users/nataliemcgee/Downloads/BedMachineGreenland-v5.nc")

# Inspect dataset contents if desired
print(ds)

bed = ds['bed'].values # Extract bed depth values
x = ds['x'].values # Extract x location (from projection)
y = ds['y'].values # Extract y location (from projection)

# Create 2D meshgrid
xx, yy = np.meshgrid(x, y)


# Define the North Polar Stereographic projection of the BedMachine data. 
# Change this if your data uses another projection
stere = pyproj.CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +datum=WGS84")
# Define the lat/lon (WGS84 geodetic) projection
geodetic = pyproj.CRS("EPSG:4326")   


# Create a transformer from stere projection to geodetic lat/lon
to_geodetic = pyproj.Transformer.from_crs(stere, geodetic, always_xy=True)

# Convert x/y to lon/lat
lon, lat = to_geodetic.transform(xx, yy)
print("Transformation Complete")

# Describe desired bounding box in lat/lon
loc_mask = (lat > 72.5) & (lat < 73.5) & (lon > -57) & (lon < -54)

# Find where mask is True
rows, cols = np.where(loc_mask)

# Get min/max indices
if len(rows) == 0 or len(cols) == 0:
    raise ValueError("No points found in the specified bounding box!")

row_min, row_max = rows.min(), rows.max()
col_min, col_max = cols.min(), cols.max()

# Slice data
lat_trim = lat[row_min:row_max+1, col_min:col_max+1]
lon_trim = lon[row_min:row_max+1, col_min:col_max+1]
bed_trim = bed[row_min:row_max+1, col_min:col_max+1]

fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.NorthPolarStereo())


#---------------------------------------------------------
# Add lat, lon points
#---------------------------------------------------------  

# Define point locations
latitudes = [
    72.91,
    72.8849,
    72.88805,
    72.8822,
    72.91523,
    72.9404,
    73.00466,
    73.017833]

longitudes = [
    -55.01,
    -55.0707,
    -55.0057,
    -54.9242,
    -55.454,
    -55.6726,
    -56.08366,
    -56.330163]

# Define section lines
vertex_lons = [-56.546, -55.308, -54.48]

vertex_lats = [73.055, 72.903, 72.89]


#---------------------------------------------------------
# Plot map
#--------------------------------------------------------- 

# Plot using PlateCarree (lon/lat) coords
pc = ax.pcolormesh(lon_trim, lat_trim, bed_trim, transform=ccrs.PlateCarree(), cmap='terrain')
ax.coastlines() # Draw coastlines
ax.set_extent([-54, -57, 72.5, 73.5], crs=ccrs.PlateCarree()) # Set map limits in lon/lat

# Add colorbar 
fig.colorbar(pc, ax=ax, label="Bed elevation (m)")

# Plot points, section lines, etc.
ax.scatter(longitudes, latitudes, color='red', s=20, transform=ccrs.PlateCarree(), label="CTD Stations")
ax.plot(vertex_lons, vertex_lats, transform=ccrs.PlateCarree(), color = 'k', ls="dashed", marker = "o", label = "Section Location")

plt.title("Upernavik Fjord\nBedMap Bathymetry and CTD Stations")
plt.legend(loc = 'lower right')
plt.show()


