import pandas as pd
import sys
from netCDF4 import Dataset
import datetime
import numpy as np
import pdb
import re
import atcf

def get_diag_name(valid_time, prefix='diag.', suffix='.nc'):
    diag_name = prefix + valid_time.strftime("%Y-%m-%d_%H.%M.%S") + suffix
    return diag_name

def origmesh(df, initfile, diagdir, debug=False):

    # Get raw values from MPAS mesh

    # input 
    # df = pandas Dataframe version of atcf data
    # init.nc = path to file with mesh cells lat/lon (first time this is run)
    #           or a dictionary containing mesh cells lat/lon (faster)
    # diagdir = path to directory with diagnostics files.


    # The first time this is called initfile is a simple string.
    # Next time, it is a dictionary with all the needed variables.
    if isinstance(initfile, str):
        if debug:
            print("reading lat/lon from", initfile)
        init = Dataset(initfile,"r")
        lonCellrad = init.variables['lonCell'][:]
        latCellrad = init.variables['latCell'][:]
        lonCell = np.degrees(lonCellrad)
        latCell = np.degrees(latCellrad)
        lonCell[lonCell >= 180] = lonCell[lonCell >=180] - 360.
        nEdgesOnCell = init.variables['nEdgesOnCell'][:]
        cellsOnCell = init.variables['cellsOnCell'][:]
        init.close()
        initfile = {
                "initfile":initfile,
                "lonCell":lonCell,
                "latCell":latCell,
                "nEdgesOnCell":nEdgesOnCell,
                "cellsOnCell":cellsOnCell
                }
    else:
        if debug:
            print("reading lat/lon from dictionary")
        lonCell = initfile["lonCell"]
        latCell = initfile["latCell"]
        nEdgesOnCell = initfile["nEdgesOnCell"]
        cellsOnCell = initfile["cellsOnCell"]
    
    itime = 0
    for index, row in df.iterrows():
        diagfile = get_diag_name(row.valid_time, prefix='diag.', suffix='.nc')

        if debug: print("reading diagfile", diagdir+diagfile)
        nc = Dataset(diagdir+diagfile, "r")

        u10  = nc.variables['u10'][itime,:]
        v10  = nc.variables['v10'][itime,:]
        mslp = nc.variables['mslp'][itime,:]
        nc.close()

        # Extract vmax, RMW, minp, and radii of wind thresholds
        raw_vmax_kts, raw_RMW_nm, raw_minp, rad_nm = atcf.derived_winds(u10, v10, mslp, lonCell, latCell, row, debug=debug)

        # TODO: figure out how to replace the row with (possibly) multiple rows with different wind radii
        # without passing df, the entire DataFrame
        df = atcf.update_df(df, row, raw_vmax_kts, raw_RMW_nm, raw_minp, rad_nm, debug=debug)
    if debug:
        print("mpas.origmesh() pausing before return")
        pdb.set_trace()
    return df, initfile


