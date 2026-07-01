#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 11:29:19 2026

@author: nataliemcgee
"""
from netCDF4 import Dataset

ncfile = Dataset("/Users/nataliemcgee/Documents/Upernavik Data/Nutrients/nitrate_profiles.nc", 'r')

# View all info
print(ncfile)

# View variable names
print(ncfile.variables.keys())

# View metadata of one variable
print(ncfile.variables["Nitrate"])

print(ncfile["Nitrate"][3][40:50])