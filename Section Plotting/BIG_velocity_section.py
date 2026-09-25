

import numpy as np
import pandas as pd

from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import gaussian_filter

from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.colors import CenteredNorm

uvel_netcdf = "/Users/nataliemcgee/Desktop/Model Output/Uvel_202108.nc"
vvel_netcdf = "/Users/nataliemcgee/Desktop/Model Output/Vvel_202108.nc"
bathymetry_file = pd.read_csv("/Users/nataliemcgee/Documents/Upernavik Data/Bathymetry Data/across_NCS_bathymetry.csv")
grid_netcdf = "/Users/nataliemcgee/Desktop/Model Output/L2_Upernavik_grid.nc"

uvel_ds = xr.open_dataset(uvel_netcdf)
vvel_ds = xr.open_dataset(vvel_netcdf)
grid_ds = xr.open_dataset(grid_netcdf)



def nearest_indices(lon2d, lat2d, lons, lats):
    tree = cKDTree(np.column_stack([lon2d.ravel(), lat2d.ravel()]))
    dist, flat_idx = tree.query(np.column_stack([lons, lats]))
    row_idx, col_idx = np.unravel_index(flat_idx, lon2d.shape)
    return row_idx, col_idx, dist


def haversine(lon1, lat1, lon2, lat2):
    R = 6371.0  # km
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)
    a = np.sin(dphi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(dlambda/2)**2
    return 2 * R * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def cumulative_distance(lons, lats):
    seg_dist = haversine(lons[:-1], lats[:-1], lons[1:], lats[1:])
    return np.concatenate([[0], np.cumsum(seg_dist)])   # km


def compute_bearing(lons, lats):
    """Local along-track bearing (radians, clockwise from north) at each
    point, using central differences (forward/backward at endpoints)."""
    lon = np.radians(lons)
    lat = np.radians(lats)
    n = len(lon)

    lon1, lat1 = np.empty(n), np.empty(n)
    lon2, lat2 = np.empty(n), np.empty(n)

    lon1[1:-1], lat1[1:-1] = lon[:-2], lat[:-2]
    lon2[1:-1], lat2[1:-1] = lon[2:], lat[2:]

    lon1[0], lat1[0] = lon[0], lat[0]
    lon2[0], lat2[0] = lon[1], lat[1]

    lon1[-1], lat1[-1] = lon[-2], lat[-2]
    lon2[-1], lat2[-1] = lon[-1], lat[-1]

    dlon = lon2 - lon1
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    bearing = np.arctan2(x, y)   # radians, clockwise from north
    return bearing


def plot_perpendicular_section(uvel_ds, vvel_ds, grid_ds, lons, lats):
    lon2d = uvel_ds["longitude"].values
    lat2d = uvel_ds["latitude"].values

    row_idx, col_idx, snap_dist = nearest_indices(lon2d, lat2d, lons, lats)
    print("max snap distance (deg):", snap_dist.max())

    # extract U and V sections, average over iterations
    u_section = uvel_ds["Uvel"].isel(
        rows=("distance", row_idx), cols=("distance", col_idx)
    ).mean(dim="iterations", skipna=True)

    v_section = vvel_ds["Vvel"].isel(
        rows=("distance", row_idx), cols=("distance", col_idx)
    ).mean(dim="iterations", skipna=True)

    # sample AngleCS/AngleSN at the same points
    angleCS = grid_ds["AngleCS"].values[row_idx, col_idx]   # shape (distance,)
    angleSN = grid_ds["AngleSN"].values[row_idx, col_idx]

    angleCS = xr.DataArray(angleCS, dims="distance")
    angleSN = xr.DataArray(angleSN, dims="distance")

    # rotate grid-relative (U, V) -> true east/north
    u_east = u_section * angleCS - v_section * angleSN
    v_north = u_section * angleSN + v_section * angleCS

    # local along-track bearing at each point (radians, clockwise from north)
    theta = compute_bearing(lons, lats)

    cos_t = xr.DataArray(np.cos(theta), dims="distance")
    sin_t = xr.DataArray(np.sin(theta), dims="distance")

    # positive v_perp = flow directed to the right of the along-track
    # direction (facing from point 0 toward the last point)  
    # perpendicular component using TRUE east/north velocities
    v_perp = u_east * cos_t - v_north * sin_t

    # real along-track distance (km), for the x-axis
    dist_km = cumulative_distance(lons, lats)

    fig, ax = plt.subplots(figsize=(20, 5))
    
    depths_arr = uvel_ds["depths"].values
    
    pcm = ax.pcolormesh(
    dist_km, depths_arr, v_perp.values,
    shading="gouraud", cmap="RdBu_r",
    norm=CenteredNorm(vcenter=0.0, halfrange=0.15),
    )

 
    section_distances = bathymetry_file['Distance_km']  
    bed = -bathymetry_file['Bed_Elevation_m'] 
    
    ax.plot(section_distances, bed, color = "gray")

    ax.fill_between(section_distances, bed, 800, color="gray", zorder=3)
    
    ax.set_xlabel("Distance along transect (km)")
    ax.set_ylabel("Depth (m)")
    ax.set_ylim(-50, 800)
    ax.set_xlim(dist_km.min(), dist_km.max())
    ax.invert_yaxis()
    ax.invert_xaxis()
    fig.colorbar(pcm, ax=ax, label="Perpendicular velocity (m/s)")
    plt.show()

    return v_perp, dist_km


# --- construct the line ---
lons1 = np.linspace(-54.57, -54.39, 200)
lats1 = np.linspace(72.845, 72.962, 200)

lons2 = np.linspace(-54.39, -54.57, 200)
lats2 = np.linspace(72.962, 72.988, 200)

# drop the first point of segment 2 since it duplicates the last point of segment 1
lons = np.concatenate((lons1, lons2[1:]))
lats = np.concatenate((lats1, lats2[1:]))

v_perp_section, dist_km = plot_perpendicular_section(uvel_ds, vvel_ds, grid_ds, lons, lats)