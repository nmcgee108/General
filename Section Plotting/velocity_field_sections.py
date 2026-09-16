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

uvel_netcdf = "/Users/nataliemcgee/Desktop/Model Output/Uvel_202107.nc"
vvel_netcdf = "/Users/nataliemcgee/Desktop/Model Output/Vvel_202107.nc"

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
    fig, ax = plt.subplots()
    pcm = ax.pcolormesh(np.arange(len(lons)), ds["depths"], section.values,
                         shading="auto", cmap="RdBu_r", norm=CenteredNorm(vcenter=0.0))
    
    ax.set_xlabel("Point index along line")
    ax.set_ylabel("Depth (m)")
    ax.set_ylim(0, 750)
    ax.set_xlim(10, 200)
    ax.invert_yaxis()
    fig.colorbar(pcm, ax=ax, label=f"{var} (m/s)")
    plt.show()

    return section

# lons = np.linspace(-54.43, -54.40, 100)
# lats = np.linspace(72.926, 72.955, 100) glacier C

# lons = np.linspace(-54.53, -54.427, 200)
# lats = np.linspace(72.98, 72.967, 200) # GLACIER N

lons = np.linspace(-54.55, -54.49, 200)
lats = np.linspace(72.850, 72.89, 200) # GLACIER S

u_section = plot_section(uvel_ds, "Uvel", lons, lats)
v_section = plot_section(vvel_ds, "Vvel", lons, lats)


