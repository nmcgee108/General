#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 09:43:58 2025

@author: nataliemcgee
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

filename = "/Users/nataliemcgee/Documents/GitHub/Upernavik-Project/Upernavik Data/ice_discharge/combined_mar_racmo_discharge1.nc"
discharge_ds = xr.open_dataset(filename)

fjords=[42704, 42733, 42968, 43247, 43649]

print(fjords)

fig, axes = plt.subplots(2, 5, sharey= True, figsize = (15, 9))

date_ranges=[["2010-01-01", "2026-01-01"],["1958-01-01", "2000-01-01"]]

for j in range(2):
    valid_dates = (discharge_ds.Date < np.datetime64(date_ranges[j][1])) & (discharge_ds.Date > np.datetime64(date_ranges[j][0]))
    for i in range(5):
        discharge_subset_racmo = discharge_ds.Discharge.sel(Fjord=fjords[i], Model="RACMO")
        discharge_subset_mar = discharge_ds.Discharge.sel(Fjord=fjords[i], Model="MAR")
        discharge_subset_racmo = discharge_subset_racmo.where(valid_dates, drop=True)
        discharge_subset_mar = discharge_subset_mar.where(valid_dates, drop=True)
        
        discharge_subset_racmo['dayofyear'] = discharge_subset_racmo.Date.dt.dayofyear
        discharge_subset_mar['dayofyear'] = discharge_subset_mar.Date.dt.dayofyear
        
        # Group by day-of-year and calculate mean and std
        daily_mean_RACMO = discharge_subset_racmo.groupby('dayofyear').mean('Date')
        daily_std_RACMO  = discharge_subset_racmo.groupby('dayofyear').std('Date')
        
        daily_mean_MAR = discharge_subset_mar.groupby('dayofyear').mean('Date')
        daily_std_MAR  = discharge_subset_mar.groupby('dayofyear').std('Date')
        
        axes[j,i].set_ylim(-5, 550)
        axes[j,i].plot(daily_mean_RACMO['dayofyear'], daily_mean_RACMO, label='RACMO mean')
        axes[j,i].fill_between(
            daily_mean_RACMO['dayofyear'],
            daily_mean_RACMO - daily_std_RACMO,
            daily_mean_RACMO + daily_std_RACMO,
            alpha=0.3, label='±1 Std Dev')
        
        axes[j,i].plot(daily_mean_MAR['dayofyear'], daily_mean_MAR, label='MAR mean')
        axes[j,i].fill_between(
            daily_mean_MAR['dayofyear'],
            daily_mean_MAR - daily_std_MAR,
            daily_mean_MAR + daily_std_MAR,
            alpha=0.3, label='±1 Std Dev')
        axes[0,0].legend()
        
        #axes[0,i].text(1,300,f"Depth: {int(discharge_subset_racmo["Depth"])}m")
    
        axes[j,i].set_xlabel("Day of Year")
        
        # Filter for 2024 only
        highlight_dates = (discharge_ds.Date >= np.datetime64("2024-01-01")) & \
                          (discharge_ds.Date < np.datetime64("2025-01-01"))
        
        # Subset for this fjord
        discharge_highlight = discharge_ds.Discharge.sel(Fjord=fjords[i])
        discharge_highlight = discharge_highlight.where(highlight_dates, drop=True)
        
        # Add day-of-year coordinate
        discharge_highlight['dayofyear'] = discharge_highlight.Date.dt.dayofyear
        
        # Plot RACMO
        axes[0, i].plot(
            discharge_highlight.sel(Model="RACMO")['dayofyear'],
            discharge_highlight.sel(Model="RACMO"),
            color="k",
            linewidth = 0.5,
            label="2013 RACMO" if i == 0 else None)
        
        # Plot MAR
        # #axes[i, 1].plot(
        #     discharge_highlight.sel(Model="MAR")['dayofyear'],
        #     discharge_highlight.sel(Model="MAR"),
        #     #linestyle="dashed",
        #     color="C1",
        #     linewidth = 1,
        #     label="2024 MAR" if i == 0 else None)
        
        # Scatter Cruise avgs
        highlight_dates = (discharge_ds.Date < np.datetime64("2024-08-02")) & (discharge_ds.Date > np.datetime64("2024-07-25"))
        #highlight_dates = (discharge_ds.Date > np.datetime64("2013-08-24")) & (discharge_ds.Date < np.datetime64("2013-09-02"))
        discharge_highlight = discharge_ds.Discharge.sel(Fjord = fjords[i])
        discharge_highlight = discharge_highlight.where(highlight_dates, drop=True)
        racmo_avg_highlight = float(discharge_highlight.mean(dim="Date")[0])
        mar_avg_highlight = float(discharge_highlight.mean(dim="Date")[1])
        
        #213 for 7/30
        axes[0, i].scatter(213, racmo_avg_highlight, color = "mediumblue", zorder = 3)
        axes[0, i].scatter(213, mar_avg_highlight, color = "orangered", zorder = 3)
        
        axes[1, i].text(1,300,f"Q_MAR: {mar_avg_highlight:.2f}")
        axes[1, i].text(1,340,f"Q_RACMO: {racmo_avg_highlight:.2f}")
        
        print()
        
        legend_elements = [Line2D([0], [0], color='C0', label='Mean RACMO +/-  1 SD'),
                           Line2D([0], [0], color='C1', label='Mean MAR +/-  1 SD'),
                           Line2D([0], [0], color='k', label='Cruise Year (RACMO)'),
                           Line2D([0], [0], marker='o', color='mediumblue', label='Cruise Week (RACMO)', linestyle='None', markersize=6),
                           Line2D([0], [0], marker='o', color='orangered', label='Cruise Week (MAR)', linestyle='None', markersize=6)]
        
        axes[0, 0].legend(handles=legend_elements, loc='best')  # or use loc='upper right', etc.

        
        
    axes[0,0].set_title("1) Unnamed", fontweight = 'bold')
    axes[0,1].set_title("2) Upernavik N", fontweight = 'bold')
    axes[0,2].set_title("3) Upernavik C", fontweight = 'bold')
    axes[0,3].set_title("4) Upernavik S", fontweight = 'bold')
    axes[0,4].set_title("5) Upernavik SS", fontweight = 'bold')
    
    axes[0,0].set_ylabel(r"$\mathbf{2010-2025}$"+"\n"+ r"Subglacial Discharge [m$^3$/s]")
    axes[1,0].set_ylabel(r"$\mathbf{1958-2000}$"+"\n"+ r"Subglacial Discharge [m$^3$/s]")
    
    

plt.tight_layout()
plt.show()

fig, axes = plt.subplots(1, figsize = (40, 5))
axes.plot(discharge_ds.Discharge.sel(Fjord=fjords[1], Model="RACMO"), label='RACMO mean')
axes.plot(discharge_ds.Discharge.sel(Fjord=fjords[2], Model="RACMO"), label='RACMO mean')
axes.plot(discharge_ds.Discharge.sel(Fjord=fjords[3], Model="RACMO"), label='RACMO mean')
axes.set_xlim(0, 5000)
