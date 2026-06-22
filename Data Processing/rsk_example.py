#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 10:44:51 2026

@author: nataliemcgee
"""

from pyrsktools import RSK


# Instantiate an RSK class object, passing the path to an RSK file

rsk_file = "/Users/nataliemcgee/Downloads/206288_20260312_1232 Deep OG March 12, 2026.rsk"

metadata = {}

with RSK(rsk_file) as rsk:

    # # Examine what attributes are populated
    #print(rsk)
    
    # # Retrieve metadata
    #print(rsk.instrument)
    #print(rsk.deployment)
    #print(rsk.channels)
    #print(rsk.schedule)
    #print(rsk.instrumentSensors)
    #print(rsk.calibrations)
    #print(rsk.epoch)
    #print(rsk.logs)
    

    # Where to find the important metadata in the rsk file: 
        
    rsk.readdata()
    data = rsk.dataArrays[0]
    print(rsk.dataArrays[0])
    
    # -------------------------------------------------------------------------
    # Dimensions
    # -------------------------------------------------------------------------
        
    # number of depth bins
    # number of profiles
    
    # -------------------------------------------------------------------------
    # Fill coordinate data
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # Fill science variable data
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # Fill or compute useful global attributes
    # -------------------------------------------------------------------------

    # title
    # summary
    # processing level
    # vessel/platform name
    
    # what instrument was used
    metadata["instrument"] = (
        getattr(rsk.instrument, "model", None)
        if getattr(rsk, "instrument", None)
        else None
    )
    
    # instrument serial number
    metadata["serial_number"] = (
        getattr(rsk.instrument, "serialID", None)
        if getattr(rsk, "instrument", None)
        else None
    )
    
    # creator name
    # creator email
    # creator institution
    # contact name
    # contact_email

    
    # geospatial lat min
    # geospatial lat max
    # geospatial lon min
    # geospatial lon max
    

    # time coverage start
    metadata["time_coverage_start"] = (
        data["timestamp"][0]
        .astype("datetime64[ms]")
        .item()
        .isoformat()
    )

    # time coverage end  
    metadata["time_coverage_end"] = (
        data["timestamp"][-1]
        .astype("datetime64[ms]")
        .item()
        .isoformat()
    )

    
print('done')
print(metadata)


















