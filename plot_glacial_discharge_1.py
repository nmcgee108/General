#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 15:29:27 2026

@author: nataliemcgee
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


plt.rcParams['font.size'] = 14

filename = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/ice_discharge/combined_mar_racmo_discharge1.nc"
discharge_ds = xr.open_dataset(filename)

fjords = [42733, 42968, 43247]

fig, axes = plt.subplots(2, 1, figsize=(9, 14))

date_ranges = [
    ["2013-01-01", "2013-12-31"],
    ["2015-01-01", "2015-12-31"],
    ["2016-01-01", "2016-12-31"],
    ["2019-01-01", "2019-12-31"],
    ["2024-01-01", "2024-12-31"]]

highlight_dates = [
    "2013-09-02", "2015-09-17", "2016-09-17",
    "2019-08-30", "2024-07-31"]

colors = plt.cm.managua(np.linspace(0, 1, 5))

for j in range(5):

    valid_dates = (
        (discharge_ds.Date <= np.datetime64(date_ranges[j][1])) &
        (discharge_ds.Date >= np.datetime64(date_ranges[j][0])))

    highlight_doy = np.datetime64(highlight_dates[j]).astype('datetime64[D]').astype(int) - \
                    np.datetime64(f"{highlight_dates[j][:4]}-01-01").astype('datetime64[D]').astype(int) + 1

    for i in range(2):

        discharge_subset = (
            discharge_ds.Discharge
            .sel(Fjord=fjords[i], Model="RACMO")
            .where(valid_dates, drop=True))

        # 7-day rolling stats
        weekly_mean = discharge_subset.rolling(Date=7, center=True, min_periods=1).mean()
        weekly_std = discharge_subset.rolling(Date=7, center=True, min_periods=1).std()

        x = weekly_mean.Date.dt.dayofyear

        upper = weekly_mean + weekly_std
        lower = weekly_mean - weekly_std

        axes[i].set_ylim(-5, 500)
        
        if j==4:
            axes[i].plot(
                x,
                weekly_mean,
                linewidth=3,
                color="darkblue",
                label=date_ranges[j][0][:4])

            axes[i].fill_between(
                x,
                lower,
                upper,
                color="darkblue",
                alpha=0.25)
        
        else: pass
            # axes[i].plot(
            #     x,
            #     weekly_mean,
            #     linewidth=2,
            #     color=colors[j],
            #     label=date_ranges[j][0][:4])
    
            # axes[i].fill_between(
            #     x,
            #     lower,
            #     upper,
            #     color=colors[j],
            #     alpha=0.25)

        # ---- Highlight cruise moment ----
        discharge_highlight = (
            discharge_ds.Discharge
            .sel(Fjord=fjords[i], Model="RACMO")
            .sel(Date=np.datetime64(highlight_dates[j])))

        racmo_avg_highlight = float(discharge_highlight)
        
        if j==4:
            axes[i].scatter(
                highlight_doy,
                racmo_avg_highlight,
                color="darkblue",
                edgecolors='black',
                s=180,
                zorder=3)
            
        else: pass
            # axes[i].scatter(
            #     highlight_doy,
            #     racmo_avg_highlight,
            #     color=colors[j],
            #     edgecolors='black',
            #     s=100,
            #     zorder=3)


        axes[i].set_xlabel("Day of Year")
        axes[i].set_ylabel(r"Subglacial Discharge [m$^3$/s]")
        
month_starts = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

for ax in axes:
    ax.set_xticks(month_starts)
    ax.set_xticklabels(month_labels)
    ax.set_xlim(125, 275)


axes[0].set_title("Upernavik N", fontweight='bold')
axes[1].set_title("Upernavik C", fontweight='bold')

#axes[1].legend(loc="upper right")

axes[0]

plt.tight_layout()
plt.show()
