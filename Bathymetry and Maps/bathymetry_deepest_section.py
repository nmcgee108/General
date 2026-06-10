#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find the deepest (maximin) route through a fjord using gridded bathymetry
based on BedMachineGreenland-v5 and plot/save a section.

Created on Mon Jun 16 17:16:03 2025
@author: nataliemcgee
Note: Modified by ChatGPT for 'deepest route' (maximizes the minimum depth along path)
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pyproj
from pyproj import Geod
import heapq

# ---------- Load data ----------
ds = xr.open_dataset("/Users/nataliemcgee/Documents/Upernavik Project/Upernavik Data/BedMachineGreenland-v5.nc")
bed = ds['bed'].values        # bed elevation in m (negative below sea level usually)
x = ds['x'].values            # projected x
y = ds['y'].values            # projected y

# Convert bed elevation to positive water depth (set negative for missing)
# If bed values are negative below sea level, depth = -bed; if bed positive (uncommon) we handle generically
depth_grid = -bed  # depth (m) positive = deeper
# mask invalid cells
depth_grid = np.where(np.isfinite(depth_grid), depth_grid, np.nan)

# Create 2D meshgrid (useful for plotting)
xx, yy = np.meshgrid(x, y)

# ---------- Projections ----------
stere = pyproj.CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +datum=WGS84")
geodetic = pyproj.CRS("EPSG:4326")  # WGS84 lat/lon
to_geodetic = pyproj.Transformer.from_crs(stere, geodetic, always_xy=True)
to_projected = pyproj.Transformer.from_crs(geodetic, stere, always_xy=True)
geod = Geod(ellps='WGS84')

# ---------- Define endpoints  ----------
start_lat, start_lon = 73.2313, -57.8811
end_lat, end_lon = 72.9339, -54.5672

# Convert start/end to projected x/y
sx, sy = to_projected.transform(start_lon, start_lat)
ex, ey = to_projected.transform(end_lon, end_lat)

# ---------- helper: find nearest grid indices for a projected point ----------
def nearest_grid_index(px, py, x_arr, y_arr):
    ix = np.argmin(np.abs(x_arr - px)) # find the index of the nearest x in the array
    iy = np.argmin(np.abs(y_arr - py)) # find the index of the nearest y in the array
    return iy, ix  # note: arrays are indexed as [iy, ix] (y first)

start_iy, start_ix = nearest_grid_index(sx, sy, x, y) # grid indeces for start point
end_iy, end_ix = nearest_grid_index(ex, ey, x, y) # grid indeces for end point

# Sanity print
print("Start grid index (y,x):", (int(start_iy), int(start_ix)), "End grid index (y,x):", (int(end_iy), int(end_ix)))
print("Start depth (m):", depth_grid[start_iy, start_ix], "End depth (m):", depth_grid[end_iy, end_ix])

# ---------- Maximin path algorithm (Dijkstra-like using a max-heap)
# We find a path from start to end that maximizes the minimum depth encountered along the path.
ny, nx = depth_grid.shape
print("Shape of Depth Grid (y, x):", (ny, nx))

# neighbor offsets (8-connectivity)
neigh_offsets = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]

# best_score stores the best achievable minimum depth to reach that cell
best_score = np.full(depth_grid.shape, -np.inf) # fill grid with - infinities
prev = np.full(depth_grid.shape + (2,), -1, dtype=int)  # create a 3d grid. store predecessor indices (two values) in each point
                                                        # (initialize filled with -1 for everything)

# initialize
best_score[start_iy, start_ix] = depth_grid[start_iy, start_ix]
# max-heap: we use negative value to use Python's min-heap as max-heap
heap = [(-best_score[start_iy, start_ix], start_iy, start_ix)]

visited = np.zeros_like(depth_grid, dtype=bool)

# The heap finds the absolute deepest "clearance" depth. 
while heap: 
    negscore, iy, ix = heapq.heappop(heap) # pull up one point to analyze and its location. heappop prioritizes the smallest "negscore". 
                                           # this will pull the path with currently the deepest (most negative) sill depth/pinch point
    score = -negscore      # score is the current shallowest depth that must be traversed to get to this point
    if visited[iy, ix]:    # skip if we've already looked at this point
        continue
    visited[iy, ix] = True
    # early stop if we reached end
    if (iy == end_iy) and (ix == end_ix):
        break
    # explore neighbors
    for dy, dx in neigh_offsets:    
        nyi, nxi = iy + dy, ix + dx
        if nyi < 0 or nxi < 0 or nyi >= ny or nxi >= nx: # quit if we hit edges of the plot
            continue
        if not np.isfinite(depth_grid[nyi, nxi]): # quit if NaN
            continue
        # candidate min depth along path to neighbor
        cand = min(score, depth_grid[nyi, nxi])  # cand = whichever is shallower, the neighbor or the current shallowest (score)
        if cand > best_score[nyi, nxi]:          # if this shallow point is deeper than the current shallowest point its a better option
            best_score[nyi, nxi] = cand          
            prev[nyi, nxi] = [iy, ix]            # update record matrix (prev) to remember how we got here
            heapq.heappush(heap, (-cand, nyi, nxi)) # store this point in the heap 

# check if we reached the end
if not np.isfinite(best_score[end_iy, end_ix]) or best_score[end_iy, end_ix] == -np.inf:
    raise RuntimeError("No path found between points with valid bathymetry.")

print("Maximin depth along best path (m):", best_score[end_iy, end_ix])

# Reconstruct path
path_indices = []
iy, ix = end_iy, end_ix     # work backwards from end
while (iy != -1) and (ix != -1):
    path_indices.append((iy, ix))  # add the point to the list
    p = prev[iy, ix]
    iy, ix = int(p[0]), int(p[1])  # for each point look in the record matrix to see how we got there
path_indices = path_indices[::-1]  # reorder start->end

# Convert path indices to projected coordinates and lat/lon
path_x = np.array([x[ix] for (iy, ix) in path_indices])
path_y = np.array([y[iy] for (iy, ix) in path_indices])
path_lon, path_lat = to_geodetic.transform(path_x, path_y)

# Interpolate depth along the path (higher resolution interpolation along the polyline)
bed_interp = RegularGridInterpolator((y, x), bed, bounds_error=False, fill_value=np.nan)
depth_along_path = -bed_interp(np.column_stack((path_y, path_x)))  # convert to positive depth

# Compute cumulative distance along path (km)
distances = [0.0]
for i in range(1, len(path_lon)):
    _, _, d = geod.inv(path_lon[i-1], path_lat[i-1], path_lon[i], path_lat[i])
    distances.append(distances[-1] + d / 1000.0)
distances = np.array(distances)

# Save CSV: Distance_km,Depth_m,Lat,Lon,proj_x,proj_y
out_arr = np.column_stack((distances, depth_along_path, path_lat, path_lon, path_x, path_y))
np.savetxt('deepest_route_section.csv', out_arr, delimiter=',',
           header='Distance_km,Depth_m,Latitude,Longitude,Proj_X,Proj_Y', comments='')
print("deepest route csv saved: deepest_route_section.csv")

# ---------- Plotting ----------
fig, axs = plt.subplots(1,2, figsize=(14,6))

# Plan view: bed (cropped) with path overlay
ax = axs[0]

# Define crop bounds with a margin (in projected coords, meters)
margin = 20000  # 20 km buffer
xmin, xmax = path_x.min() - margin, path_x.max() + margin
ymin, ymax = path_y.min() - margin, path_y.max() + margin

# Select only the indices within the bounding box
ix = np.where((x >= xmin) & (x <= xmax))[0]
iy = np.where((y >= ymin) & (y <= ymax))[0]

# Make cropped grids
x_sub = x[ix]
y_sub = y[iy]
depth_sub = depth_grid[np.ix_(iy, ix)]
xx_sub, yy_sub = np.meshgrid(x_sub, y_sub)

# Plot cropped bathymetry
c = ax.pcolormesh(xx_sub, yy_sub, depth_sub, cmap = "terrain", shading='auto')

# Overlay path
ax.plot(path_x, path_y, '-r', linewidth=0.5, label='Deepest path')
ax.scatter([sx, ex], [sy, ey], c='y', marker='*', zorder=10, label='start/end')

ax.set_title("Plan view (cropped) with deepest path")
ax.set_aspect('equal')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.legend(loc="lower right")
fig.colorbar(c, ax=ax, label='Depth (m)')


# Cross-section (distance vs depth)
ax2 = axs[1]
ax2.plot(distances, depth_along_path, color='steelblue')
ax2.fill_between(distances, depth_along_path, 1000)
ax2.set_ylim(0,1000)
ax2.set_xlabel("Distance along path (km)")
ax2.set_ylabel("Depth (m)")
ax2.invert_yaxis()  # if you like depth downwards (optional)
ax2.set_title("Bathymetric section along deepest path")

plt.tight_layout()
plt.show()

# Also optionally save a KML or a small ascii summary
print("Path length (km):", distances[-1])
print("Minimum depth along path (m):", np.nanmin(depth_along_path))
print("Maximum depth along path (m):", np.nanmax(depth_along_path))
