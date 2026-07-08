#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 11:29:19 2026

@author: nataliemcgee
"""
from netCDF4 import Dataset

ncfile = Dataset("/Users/nataliemcgee/Documents/Upernavik Data/Morven CTD Data 2013-2019/2015_profiles.nc", 'r')

# View all info
print(ncfile)

# View variable names
print(ncfile.variables.keys())

# View metadata of one variable
print(ncfile.variables["Conservative_Temperature"])

print(ncfile["Nitrate"][3][40:50])