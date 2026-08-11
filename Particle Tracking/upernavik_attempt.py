#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 10:22:24 2026

@author: nataliemcgee
"""

import math
from datetime import timedelta
from operator import attrgetter

import matplotlib.pyplot as plt
import numpy as np
import trajan as ta
import xarray as xr
from IPython.display import HTML
from matplotlib.animation import FuncAnimation

import warnings
import parcels

warnings.filterwarnings("ignore", category=UserWarning, module="parcels")


filenames = {
    "U": "/Users/nataliemcgee/Desktop/Model Output/Uvel_202106.nc",
    "V": "/Users/nataliemcgee/Desktop/Model Output/Vvel_202106.nc",
}

variables = {
    "U": "Uvel",
    "V": "Vvel",
}
dimensions = {"lat": "latitude", "lon": "longitude", "time": "iterations"}

fieldset = parcels.FieldSet.from_netcdf(filenames, variables, dimensions)


print(fieldset)

pset = parcels.ParticleSet.from_line(
    fieldset=fieldset,
    pclass=parcels.JITParticle,
    size=5,  # releasing 5 particles
    start=(-55.5, 73),  # releasing on a line: the start longitude and latitude
    finish=(-56.5, 73),  # releasing on a line: the end longitude and latitude
)

print(pset)

output_file = pset.ParticleFile(
    name="Upernavik_ex.zarr", outputdt=timedelta(hours=6)
)
pset.execute(
    parcels.AdvectionRK4,
    runtime=timedelta(days=0.001),
    dt=timedelta(minutes=5),
    output_file=output_file,
)

ds = xr.open_zarr("Upernavik_ex.zarr")
ds.traj.plot(margin=2)
plt.show()










