#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 13:48:32 2026

@author: nataliemcgee
"""
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.colors import CenteredNorm

uvel_netcdf = "/Users/nataliemcgee/Desktop/Model Output/Uvel_202108.nc"
vvel_netcdf = "/Users/nataliemcgee/Desktop/Model Output/Vvel_202108.nc"

uvel_ds = xr.open_dataset(uvel_netcdf)
vvel_ds = xr.open_dataset(vvel_netcdf)

def plot_section(ds, var, lons, lats):
    lon2d = ds["longitude"].values
    lat2d = ds["latitude"].values

    # nearest grid cell for each point on the line
    tree = cKDTree(np.column_stack([lon2d.ravel(), lat2d.ravel()]))
    _, flat_idx = tree.query(np.column_stack([lons, lats]))
    row_idx, col_idx = np.unravel_index(flat_idx, lon2d.shape)
    
    

    # pull out the transect (dims: depths, distance)
    section = ds[var].isel(
        rows=("distance", row_idx),
        cols=("distance", col_idx)
    ).mean(dim="iterations")

    # plot using point index as x-axis (skip real distance calc for simplicity)
    fig, ax = plt.subplots(figsize = (10, 5))
    pcm = ax.pcolormesh(np.arange(len(lons)), ds["depths"], section.values,
                         shading="auto", cmap="RdBu_r", norm=CenteredNorm(vcenter=0.0))
    
    ax.set_xlabel("Point index along line")
    ax.set_ylabel("Depth (m)")
    ax.set_ylim(0, 800)
    ax.set_xlim(0, 400)
    ax.invert_yaxis()
    fig.colorbar(pcm, ax=ax, label=f"{var} (m/s)")
    plt.show()

    return section

# lons = np.linspace(-54.43, -54.40, 100)
# lats = np.linspace(72.926, 72.955, 100) glacier C

# lons = np.linspace(-54.53, -54.427, 200)
# lats = np.linspace(72.98, 72.967, 200) # GLACIER N

# lons = np.linspace(-54.55, -54.49, 200)
# lats = np.linspace(72.850, 72.89, 200) # GLACIER S

lons1=np.linspace(-54.57, -54.39, 200)
lons2=np.linspace(-54.39, -54.57, 200)

lats1=np.linspace(72.845, 72.962, 200)
lats2=np.linspace(72.962, 72.988, 200)

lons = np.concat((lons1, lons2))
lats = np.concat((lats1, lats2))


u_section = plot_section(uvel_ds, "Uvel", lons, lats)
v_section = plot_section(vvel_ds, "Vvel", lons, lats)


