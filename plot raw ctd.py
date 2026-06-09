#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 11:30:54 2025

@author: nataliemcgee
"""
from read_cnv_header import read_cnv_header
from parse_cnv import parse_cnv
import matplotlib.pyplot as plt

filename = '/Users/nataliemcgee/Documents/Upernavik Project/Upernavik Data/cnv_data/2024-07-31T115713_SBE0251108_2.cnv'

#Parse header for variable names and other information
cast_num, latitude, longitude, variables = read_cnv_header(filename)

#Extract data from cnv file
#data_cols is a dictionary with key=variable names
#data_array is an array with each column = one variable (corresponds with index)
data_cols, data_array= parse_cnv(filename, variables)

#Grab individual variable data from data_cols as needed
#Use variable name (from print statement) as key: i.e. ['timeS: Time, Elapsed [seconds]']

pres = data_cols['prdM: Pressure, Strain Gauge [db]']
oxy = data_cols['sbeox0ML/L: Oxygen, SBE 43 [ml/l]']
fluor = data_cols['flECO-AFL: Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]']
turb2 = data_cols['turbWETntu1: Turbidity, WET Labs ECO, 2 [NTU]']
cond = data_cols['c0mS/cm: Conductivity [mS/cm]']

plt.figure(figsize = (4, 8))
plt.plot(oxy, -pres)
plt.xlim(4, 11)
plt.ylim(-300, -200)
plt.xlabel('Oxygen')
plt.ylabel('(Negative) Pressure')
plt.title("Cast 2")
