#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 12:27:38 2025

@author: nataliemcgee
"""

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pyproj
import numpy as np
from scipy.signal import savgol_filter

plt.rcParams['font.size'] = 18

# Load the NetCDF dataset
ds = xr.open_dataset("/Users/nataliemcgee/Documents/Upernavik Project/Upernavik Data/BedMachineGreenland-v5.nc")

# Inspect dataset contents
print(ds)

bed = ds['bed'].values
x = ds['x'].values
y = ds['y'].values


# Create 2D meshgrid
xx, yy = np.meshgrid(x, y)

# Define polar stereographic projection used by BedMachine Greenland
#proj = pyproj.Proj('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +datum=WGS84')
proj = pyproj.Proj('EPSG:3413')  # NSIDC Sea Ice Polar Stereographic North

# Convert x/y to lon/lat
lon, lat = proj(xx, yy, inverse=True)
print("done")

# find bounding box in lat/lon
loc_mask = (lat > 72.6) & (lat < 75.5) & (lon > -63.5) & (lon < -55)

# find where mask is True
rows, cols = np.where(loc_mask)

# get min/max indices
if len(rows) == 0 or len(cols) == 0:
    raise ValueError("No points found in the specified bounding box!")

row_min, row_max = rows.min(), rows.max()
col_min, col_max = cols.min(), cols.max()

# slice data
lat_trim = lat[row_min:row_max+1, col_min:col_max+1]
lon_trim = lon[row_min:row_max+1, col_min:col_max+1]
bed_trim = bed[row_min:row_max+1, col_min:col_max+1]



# ## TEMPORARY: SAVE LOCAL DATA
# np.savez(
#     "melvillebay_bedmachine_subset1.npz",
#     bed=bed_trim,
#     lat=lat_trim,
#     lon=lon_trim
# )


fig = plt.figure(figsize=(15, 6))
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-55))

# Plot using PlateCarree (lon/lat) coords
bed_trim_smooth = savgol_filter(bed_trim, window_length=7, polyorder=2, axis=0)
pc = ax.pcolormesh(lon_trim, lat_trim, bed_trim, transform=ccrs.PlateCarree(), cmap='gist_gray')
cs = ax.contour(lon_trim, lat_trim, bed_trim_smooth, 
                 levels=[-600, -400, -200],
                 colors=['mediumpurple', 'teal', 'yellow'],
                 transform=ccrs.PlateCarree())

manual_locs = [(-96890.61, -1912000), (-97324.61, -1907511.40), (-99092.62, -1900133.34)]
ax.clabel(cs, fmt="%.0f", manual=manual_locs, fontsize = 14)
#ax.clabel(cs, cs.levels, fontsize=14, inline=True, inline_spacing=10)
#[(np.float64(-96890.61492485134), np.float64(-1912285.446149565))
# (np.float64(-97324.61871355172), np.float64(-1907511.4044738603)),
# (np.float64(-98192.62629095261), np.float64(-1900133.340065953))]



#ax.coastlines() # draw coastlines


ax.set_extent([-58.5, -54, 73.4, 72.65], crs=ccrs.PlateCarree()) # Set map limits in lon/lat


# Add colorbar 
cbar = fig.colorbar(pc, ax=ax, pad=0.02)
cbar.set_label("Bed Elevation [m]")

#fjord
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

#trough
# latitudes = [
#     72.7735,
#     72.83995,
#     72.91643,
#     72.9867,
#     73.05838,
#     73.12492,
#     73.19801,
#     73.27043,
#     73.3421]

# longitudes = [
#     -58.58975,
#     -58.7702,
#     -58.92601,
#     -59.146,
#     -59.34056,
#     -59.504433,
#     -59.70404,
#     -59.90671,
#     -60.076983]


startlat, startlon = 73.24633333333334, -58.31155555555556
midlat, midlon = 72.89061111111111, -55.23407777777778
endlat, endlon = 72.89, -54.7


ax.scatter(longitudes, latitudes, color='red', s=40, transform=ccrs.PlateCarree(), zorder = 4, label="CTD Stations")
ax.scatter(-57.3294, 73.127, color='darkorange', s=120, marker = "*", transform=ccrs.PlateCarree(), zorder = 4, label="Shelf Station 9")

ax.plot([startlon, midlon], [startlat, midlat], transform=ccrs.PlateCarree(), color = 'white', ls="dashed", marker = "o")
ax.plot([midlon, endlon], [midlat, endlat], transform=ccrs.PlateCarree(), color = 'white', ls="dashed", marker = "o")

# vertex_lons = [-58.58975,-60.3796]
# vertex_lats = [72.7735, 73.4532]

plt.legend(loc = 'lower left')
plt.show()


