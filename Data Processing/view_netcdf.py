#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 11:29:19 2026

@author: nataliemcgee
"""
from netCDF4 import Dataset

ncfile = Dataset("/Users/nataliemcgee/Downloads/SANNA_2016_SBE19plus_CTD_profiles.nc", 'r')

# View all info
print(ncfile)

# View variable names
print(ncfile.variables.keys())

# View metadata of one variable
print(ncfile.variables["latitude"][0])
print(ncfile.variables["longitude"][0])


import csv

lats = ncfile.variables["latitude"][0]
lons = ncfile.variables["longitude"][0]

with open("latlon.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["latitude", "longitude"])   # header
    for lat, lon in zip(lats, lons):
        writer.writerow([lat, lon])
#         writer.writerow([lat, lon])
        
#print(ncfile["Nitrate"][3][40:50])