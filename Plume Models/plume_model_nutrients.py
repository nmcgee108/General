#!/usr/bin/env python3

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from full_plume import run_plume
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)

ctd_netcdf       = "/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc"
nutrients_file   = pd.read_csv("/Users/nataliemcgee/Documents/Upernavik Data/Nutrients/NutrientsUS2024_plotting.csv", encoding="latin-1")
nutrient_profiles = "/Users/nataliemcgee/Documents/Upernavik Data/Nutrients/nitrate_profiles.nc"

nitrate_ds = xr.open_dataset(nutrient_profiles)
ctd_ds     = xr.open_dataset(ctd_netcdf)


def run_model(cast_num, plume_depth):

    Q_discharge = 55   # subglacial discharge (m3/s)
    width       = 300  # (m)

    ctd_depth = ctd_ds["depth"].values
    ctd_sal   = ctd_ds["SAL_ABSOLUTE"][cast_num - 1].values
    ctd_temp  = ctd_ds["CONSERVATIVE_TEMP"][cast_num - 1].values

    depth   = nitrate_ds["depth"].values
    nitrate = nitrate_ds["Nitrate"][cast_num - 1].values

    # Pad CTD arrays to match nitrate depth length
    target_len   = nitrate.shape[0]
    sal_padded   = np.pad(ctd_sal,  (0, target_len - ctd_sal.shape[0]),  constant_values=np.nan)
    temp_padded  = np.pad(ctd_temp, (0, target_len - ctd_temp.shape[0]), constant_values=np.nan)

    valid = np.where(~np.isnan(ctd_sal))[0]
    ctd_maxdepth = float(ctd_depth[valid[-1]]) if len(valid) > 0 else 0.0

    # Build mask from `depth` (nitrate grid) so all arrays index consistently
    Q0         = Q_discharge / width
    depth_mask = np.where(depth < plume_depth)

    zi = -depth[depth_mask][::-1]   # deepest to shallowest, negative convention
    xi = np.zeros_like(zi)
    Ta = temp_padded[depth_mask][::-1]
    Sa = sal_padded[depth_mask][::-1]
    Na = nitrate[depth_mask][::-1]
    alpha = 0.1

    if ctd_maxdepth < plume_depth:
        print("****Warning: Chosen plume depth exceeds CTD data! T, S extrapolated at depth*****")
        zi = np.insert(zi, 0, -plume_depth)
        xi = np.insert(xi, 0, 0)
        Ta = np.insert(Ta, 0, np.nanmean(Ta[:10]))
        Sa = np.insert(Sa, 0, np.nanmean(Sa[:10]))
        Na = np.insert(Na, 0, np.nanmean(Na[:10]))

    sol = run_plume(zi, xi, Ta, Sa, Na, Q0, alpha)

    NBD             = sol["zNB"]
    nitrate_NBD     = sol["NNB"]
    temp_NBD        = sol["TNB"]
    sal_NBD        = sol["SNB"]
    amb_nitrate_NBD = nitrate[-int(sol["zNB"])]
    nitrate_anomaly = nitrate_NBD - amb_nitrate_NBD
    volume_flux     = sol["QNB"] * width
    NFA             = (nitrate_NBD - amb_nitrate_NBD) * volume_flux / 1e3

    # print(f"NBD = {NBD:.2f} [meters]")
    # print(f"Volume Flux = {volume_flux:.2f} [m3/s]")
    print(f"Nitrate @ NBD = {nitrate_NBD:.2f} [uM]")
    # print(f"Ambient Nitrate @ NBD = {amb_nitrate_NBD:.2f} [uM]")
    print(f"Nitrate Anomaly @ NBD = {nitrate_anomaly:.2f} [uM]")
    # print(f"***** NFA = {NFA:.2f} [mol/s] *****")

    return nitrate_NBD, temp_NBD, sal_NBD, NBD


# ── Section plot ──────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 8), sharey=True)
axes[0].set_ylabel("Depth [m]")
axes[0].set_xlabel(r"Nitrate [$\mu$M]")
axes[1].set_xlabel(r"CT [°C]")

depth     = nitrate_ds["depth"].values
nitrate   = nitrate_ds["Nitrate"].values
castnums  = nitrate_ds["cast"].values
sal       = ctd_ds["SAL_ABSOLUTE"].values
temp      = ctd_ds["CONSERVATIVE_TEMP"].values
ctd_depth = ctd_ds["depth"].values

for i in [1, 2, 13]:
    if castnums[i] in [1, 4, 5, 10, 11, 12]:
        continue

    zorder    = 1
    linewidth = 2
    label     = None

    if castnums[i] < 9:
        color  = "indigo"
        zorder = 2
        if castnums[i] == 3:
            label = "C. 2 and 3"

    if castnums[i] == 9:
        color     = "k"
        zorder    = 3
        linewidth = 3
        label     = "C. 9"

    if castnums[i] > 9:
        color     = "gray"
        linewidth = 1
        if castnums[i] == 14:
            label = "Trough"

    if castnums[i] == 14:
        color     = "red"
        linewidth = 2
        zorder    = 3
        label     = "C. 14"

    axes[0].plot(nitrate[i], -depth, color=color, zorder=zorder,
              linewidth=linewidth, label=label)
    
    axes[1].plot(temp[i], -ctd_depth, color=color, zorder=zorder,
              linewidth=linewidth, label=label)
    

# Plume model results — scatter at (nitrate @ NBD, NBD depth)


for cast, pdepth, marker, color, label in [(3, 175, "o", "indigo", "Using cast 3"), 
                                    (14, 175, "o", "red", "Using cast 14"),
                                    (3, 650, "s", "indigo", ""), 
                                    (14, 650, "s", "red", "")]:
    
    nitrate_NBD, temp_NBD, sal_NBD, NBD = run_model(cast, pdepth)
    
    axes[0].scatter(nitrate_NBD, NBD, marker=marker, zorder=5,
                 color = color, label = label, edgecolor="k", s=75)
    
    axes[1].scatter(temp_NBD, NBD, marker=marker, zorder=5,
                 color = color, label = label, edgecolor="k", s=75)
    
    
    
axes[0].set_xlim(7, 20)
axes[1].set_xlim(-1.5, 4)
axes[1].set_ylim(-600, 50)

axes[0].legend(loc="lower left")
plt.tight_layout()
plt.show()






