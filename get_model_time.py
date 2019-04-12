from netCDF4 import Dataset
import datetime
import pytz
import pdb
import re
import sys
import os

def valid(ncfilename, diagnostic_name):

    nc = Dataset(ncfilename,"r")
    try:
        x = nc.variables[diagnostic_name]
    except KeyError:
        print(diagnostic_name, "not found. Choices:", list(nc.variables.keys()))
        sys.exit(1)
    global_atts = nc.ncattrs()
    if 'valid_date' in global_atts:
        valid_time = nc.valid_date
        # convert unicode to datetime object
        valid_time = datetime.datetime.strptime(valid_time, '%Y-%m-%d_%H:%M:%S')
    elif 'initial_time' in x.ncattrs():
        initialization_time = datetime.datetime.strptime(x.initial_time, '%m/%d/%Y (%H:%M)')
        if x.forecast_time_units == 'hours':
            td = datetime.timedelta(0,0,0,0,0,1)
        if x.forecast_time_units == '15 minutes':
            td = datetime.timedelta(0,0,0,0,15,0)
        if x.forecast_time_units == '30 minutes':
            td = datetime.timedelta(0,0,0,0,30,0)

        forecast_lead_time = x.forecast_time * td
        valid_time = initialization_time + forecast_lead_time
    elif 'valid_time' in x.ncattrs():
        valid_time = datetime.datetime.strptime(x.valid_time, '%Y%m%d_%H%M%S')
        initialization_time = datetime.datetime.strptime(x.init_time, '%Y%m%d_%H%M%S')
    elif 'START_DATE' in global_atts:
        # Like ds300 NCAR WRF ensemble diags files
        initialization_time = datetime.datetime.strptime(nc.START_DATE, '%Y-%m-%d_%H:%M:%S')
        basename = os.path.basename(ncfilename)
        # start with "diag", then anything, then _d[0-9][0-9], then _yyyymmddhh, then ...
        m = re.search('diag.*_d\d\d_201\d{7}.*_f(\d\d\d).nc', basename)
        if m:
            fhr_str = m.group(1)
            fhr = int(fhr_str)
            forecast_lead_time = datetime.timedelta(hours=fhr)
            valid_time = initialization_time + forecast_lead_time
        else:
            print("Could not get valid time from "+basename+". Using init time instead")
            valid_time = initialization_time
    else:
        print("don't know how to get time from", ncfilename)
        print(x)
        pdb.set_trace()

    nc.close()

    # return time_zone_aware datetimes
    return pytz.utc.localize(valid_time), pytz.utc.localize(initialization_time)

def init(ncfilename, diagnostic_name):
    valid_time, init_time = valid(ncfilename, diagnostic_name)
    return init_time

