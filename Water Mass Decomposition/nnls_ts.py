#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 14:33:08 2026

@author: nataliemcgee
"""

from scipy.optimize import nnls
import xarray as xr
import numpy as np
import gsw
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 16

# ── Load data ─────────────────────────────────────────────────────────────────
ctd_netcdf = ("/Users/nataliemcgee/Documents/Upernavik Data/Padded CTD Datasets/uc_patch_dataset_padded.nc")

ctd_ds = xr.open_dataset(ctd_netcdf)

ctd_sal      = ctd_ds["SAL_ABSOLUTE"].values        # (casts, depths)
ctd_pres     = ctd_ds["PRESSURE"].values
ctd_temp     = ctd_ds["CONSERVATIVE_TEMP"].values
ctd_oxy      = ctd_ds["OXYGEN"].values
ctd_lats     = ctd_ds["LAT"].values
ctd_lons     = ctd_ds["LON"].values

sigma0        = gsw.sigma0(ctd_sal, ctd_temp)
converted_oxy = ctd_oxy * 44.660 / ((1000 + sigma0) / 1000)  # mL/L → µmol/kg

# ── End-member T / S / O values ───────────────────────────────────────────────
T_AW,  S_AW,  O_AW  =  2.85, 34.6,  233.10
T_PW,  S_PW,  O_PW  = 0, 34.1, 300 #0.55, 33.6, 270.13#-1.5,  33.2,  326.78
T_SMW, S_SMW, O_SMW = -87.0,  0.0, 1050.0
T_SGD, S_SGD, O_SGD =   0.0,  0.0,  457.0

WM_INDEX = {"AW": 0, "PW": 1, "SMW": 2, "SGD": 3}  # used in solver + plot

# ── Build & condition the OMP system ─────────────────────────────────────────
A = np.array([
    [T_AW,  T_PW,  T_SMW,  T_SGD],
    [S_AW,  S_PW,  S_SMW,  S_SGD],
    [O_AW,  O_PW,  O_SMW,  O_SGD],
    [1,     1,     1,      1    ]
], dtype=float)

Amean   = A[:-1].mean(axis=1)         # (3,) mean  over end-members, per property
Astd    = A[:-1].std(axis=1, ddof=0)  # (3,) stdev over end-members, per property

A_norm  = np.vstack(
    [(A[:-1] - Amean[:, None]) / Astd[:, None],  # normalised T, S, O rows
     np.ones((1, 4))]                              # mass-conservation row
)

Wx      = np.array([121, 97, 1, 121], dtype=float)  # Lindeman et al. 2024 weights
Aw_norm = Wx[:, None] * A_norm                       # (4, 4) weighted system

print(f"Condition number  A      : {np.linalg.cond(A):.3e}")
print(f"Condition number  Aw_norm: {np.linalg.cond(Aw_norm):.3e}")

# ── Mixing line helpers ───────────────────────────────────────────────────────
T_deepwater = 2.550228
S_deepwater = 34.697933
T_ice_eff   = -90.0

def runoff_line(s):
    return (T_deepwater / S_deepwater) * s

def melting_line(s):
    return ((T_deepwater - T_ice_eff) / S_deepwater) * s + T_ice_eff

# ── Vectorised OMP solver ─────────────────────────────────────────────────────
def compute_fractions(cast_range):
    """
    Solve NNLS for every valid observation across *all* casts in one pass.

    Returns
    -------
    fractions : ndarray (N, 4)  — AW / PW / SMW / SGD fractions [0–1]
    S_all     : ndarray (N,)
    T_all     : ndarray (N,)
    """
    # ── Collect all profiles into flat 1-D arrays ──────────────────────────
    S_all = np.concatenate([ctd_sal[c - 1]        for c in cast_range])
    T_all = np.concatenate([ctd_temp[c - 1]       for c in cast_range])
    O_all = np.concatenate([converted_oxy[c - 1]  for c in cast_range])
    n     = len(S_all)

    # ── Normalise all points at once (broadcast, no Python loop) ──────────
    obs      = np.column_stack([T_all, S_all, O_all])   # (N, 3)
    obs_norm = (obs - Amean) / Astd                       # (N, 3)  broadcast

    # Append mass-conservation column, then apply weights → (N, 4)
    dw_norm = np.hstack([obs_norm, np.ones((n, 1))]) * Wx

    # ── Only run NNLS on rows that have no NaN ─────────────────────────────
    valid     = np.all(np.isfinite(obs), axis=1)          # (N,) bool
    fractions = np.full((n, 4), np.nan)

    for i in np.flatnonzero(valid):                       # tight C-level loop
        fractions[i], _ = nnls(Aw_norm, dw_norm[i])

    return fractions, S_all, T_all


# ── Plotting ──────────────────────────────────────────────────────────────────
def plot_TS(water_mass, casts, min_conc, max_conc, colormap):
    """
    T/S diagram coloured by OMP water-mass fraction.

    Parameters
    ----------
    water_mass : "AW" | "PW" | "SMW" | "SGD"
    casts      : (first_cast, last_cast)  inclusive, 1-based
    min_conc   : colorbar minimum [%]
    max_conc   : colorbar maximum [%]
    colormap   : matplotlib colormap name
    """
    cast_range = range(casts[0], casts[1] + 1)
    wm_idx     = WM_INDEX[water_mass]

    # ── Solve OMP for every cast × depth in one call ───────────────────────
    fractions, S_all, T_all = compute_fractions(cast_range)
    frac_pct = fractions[:, wm_idx] * 100          # → percentage

    # ── Isopycnal grid (visible range only — no wasted computation) ────────
    sx = np.arange(30.0, 40, 0.1)
    ty = np.arange(-2.0,  4.05, 0.1)
    S_grid, T_grid = np.meshgrid(sx, ty)
    PDEN = gsw.rho(S_grid, T_grid, 0) - 1000

    fig, ax = plt.subplots(figsize=(10, 8))

    # Manually label contours
    contour = ax.contour(S_grid, T_grid, PDEN, levels=[23, 23.5, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28], colors='grey')
    label_positions = [(30.5, 4), (30.5, 2), (31.2, 1.5), (32, 2), (32.6, 2), (33.2, 2.1), (33.8, 2.2), (34.5, 1), (35.0, 2)]
    plt.clabel(contour, inline=True, manual=label_positions, fmt='%1.1f')

    # ── Mixing lines ───────────────────────────────────────────────────────
    sx_line = np.linspace(30, 35, 300)
    ax.plot(sx_line, runoff_line(sx_line),
            color='orange', ls='--', label='Runoff line',  zorder=2)
    ax.plot(sx_line, melting_line(sx_line),
            color='pink',   ls='--', label='Melting line', zorder=2)

    # ── Single scatter call ─────────────────
    # Pre-filter to the visible range: avoids matplotlib clipping thousands
    # of invisible points, and keeps the PathCollection lean.
    in_view = (
        (S_all >= 32.5) & (S_all <= 35) &
        (T_all >= -2) & (T_all <=  4) &
        np.isfinite(frac_pct)
    )

    sc = ax.scatter(
        S_all[in_view], T_all[in_view],
        c        = frac_pct[in_view],
        cmap     = colormap,
        vmin     = min_conc,
        vmax     = max_conc,
        s        = 5,
        zorder   = 2,
        rasterized = True,
    )

    cbar = fig.colorbar(sc, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label(f"{water_mass} fraction [%]")

    ax.set_xlabel("Absolute Salinity [g/kg]")
    ax.set_ylabel("Conservative Temperature [°C]")
    ax.set_xlim(32.5, 34.75)
    ax.set_ylim(-0.25,  3.25)
    ax.legend(loc='upper left')
    ax.set_title(f"{water_mass} Fraction", fontweight='bold')

    plt.tight_layout()
    return fig, ax


# ── Run ───────────────────────────────────────────────────────────────────────
fig, ax = plot_TS(
    water_mass = "AW",
    casts      = (2, 8),
    min_conc   = 0,
    max_conc   = 100,
    colormap   = "YlOrRd",
)
plt.show()

fig, ax = plot_TS(
    water_mass = "PW",
    casts      = (2, 9),
    min_conc   = -10,
    max_conc   = 50,
    colormap   = "BuPu",
)
plt.show()


fig, ax = plot_TS(
    water_mass = "SMW",
    casts      = (2, 8),
    min_conc   = 0,
    max_conc   = 1,
    colormap   = "PuBu",
)
plt.show()


fig, ax = plot_TS(
    water_mass = "SGD",
    casts      = (2, 8),
    min_conc   = 0,
    max_conc   = 1,
    colormap   = "YlGn",
)
plt.show()










