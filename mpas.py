import datetime
from fieldinfo import fieldinfo,readNCLcm
import logging
import numpy as np
import os
import re
import sys

# fieldinfo should have been imported from fieldinfo module.
# Copy fieldinfo dictionary for MPAS. Change some fnames and filenames.
fieldinfo['precip']['fname'] = ['rainnc']
fieldinfo['precip-24hr']['fname'] = ['rainnc']
fieldinfo['precip-48hr']['fname'] = ['rainnc']
fieldinfo['precipacc']['fname'] = ['rainnc']
fieldinfo['precipacc']['filename'] = 'diag'
fieldinfo['sbcape']['fname'] = ['sbcape']
fieldinfo['sbcape']['filename'] = 'diag'
fieldinfo['mlcape']['fname'] = ['mlcape']
fieldinfo['mlcape']['filename'] = 'diag'
fieldinfo['mucape']['fname'] = ['cape']
fieldinfo['mucape']['filename'] = 'diag'
fieldinfo['sbcinh']['fname'] = ['sbcin']
fieldinfo['sbcinh']['filename'] = 'diag'
fieldinfo['mlcinh']['fname'] = ['mlcin']
fieldinfo['mlcinh']['filename'] = 'diag'
fieldinfo['pwat']['fname'] = ['precipw']
fieldinfo['pwat']['filename'] = 'diag'
fieldinfo['mslp']['fname'] = ['mslp']
fieldinfo['mslp']['filename'] = 'diag'
fieldinfo['td2']['fname'] = ['surface_dewpoint']
fieldinfo['td2depart']['fname'] = ['surface_dewpoint']
fieldinfo['thetae']['fname'] = ['t2m', 'q2', 'surface_pressure']
fieldinfo['rh2m']['fname'] = ['t2m', 'surface_pressure', 'q2']
fieldinfo['pblh']['fname'] = ['hpbl']
fieldinfo['hmuh']['fname'] = ['updraft_helicity_max']
fieldinfo['hmuh03']['fname'] = ['updraft_helicity_max03']
fieldinfo['hmuh01']['fname'] = ['updraft_helicity_max01']
fieldinfo['rvort1']['fname'] = ['rvort1_max']
fieldinfo['hmup']['fname'] = ['w_velocity_max']
fieldinfo['hmdn']['fname'] = ['w_velocity_min']
fieldinfo['hmwind']['fname'] = ['wind_speed_level1_max']
fieldinfo['hmgrp']['fname'] = ['grpl_max']
fieldinfo['cref']['fname'] = ['refl10cm_max']
fieldinfo['cref']['filename'] = 'diag'
fieldinfo['ref1km']['fname'] = ['refl10cm_1km']
fieldinfo['ref1km']['filename'] = 'diag'
for ztop in ['3','1']:
    fieldinfo['srh'+ztop]['fname'] = ['srh_0_'+ztop+'km']
    fieldinfo['srh'+ztop]['filename'] =  'diag'
for ztop in ['6','1']:
    fieldinfo['shr0'+ztop+'mag']['fname'] = ['uzonal_'+ztop+'km', 'umeridional_'+ztop+'km', 'uzonal_surface', 'umeridional_surface']
    fieldinfo['shr0'+ztop+'mag']['filename'] = 'diag'
fieldinfo['zlfc']['fname'] = ['lfc']
fieldinfo['zlcl']['fname'] = ['lcl']
fieldinfo['zlcl']['filename'] = 'diag' # only zlcl filename needed to be changed from upp, not zlfc
for plev in ['200', '250','300','500','700','850','925']:
    fieldinfo['hgt'+plev]['fname'] = ['height_'+plev+'hPa']
    fieldinfo['speed'+plev]['fname'] = ['uzonal_'+plev+'hPa','umeridional_'+plev+'hPa']
    del fieldinfo['speed'+plev]['arraylevel']
    fieldinfo['temp'+plev]['fname'] = ['temperature_'+plev+'hPa']
    del fieldinfo['temp'+plev]['arraylevel']
for plev in ['500', '700', '850']:
    fieldinfo['td'+plev]['fname'] = ['dewpoint_'+plev+'hPa']
    del fieldinfo['td'+plev]['arraylevel']
    fieldinfo['vort'+plev] = {'levels' : [0,9,12,15,18,21,24,27,30,33], 'cmap': readNCLcm('prcp_1'), 'fname': ['vorticity_'+plev+'hPa'], 'filename':'diag'}
fieldinfo['vortpv']        = {'levels' : [0,9,12,15,18,21,24,27,30,33], 'cmap': readNCLcm('prcp_1'), 'fname': ['vort_pv'],               'filename':'diag'}
for plev in ['300', '500', '700', '850', '925']:
    fieldinfo['rh'+plev]['fname'] = ['relhum_'+plev+'hPa']
fieldinfo['speed10m']['fname'] = ['u10', 'v10']
fieldinfo['speed10m-tc']['fname'] = ['u10','v10']
fieldinfo['stp']['fname'] = ['sbcape','lcl','srh_0_1km','uzonal_6km','umeridional_6km','uzonal_surface','umeridional_surface']
fieldinfo['stp']['filename'] = 'diag'
fieldinfo['crefuh']['fname'] = ['refl10cm_max', 'updraft_helicity_max']
fieldinfo['crefuh']['filename'] = 'diag'
fieldinfo['wind10m']['fname'] = ['u10','v10']
fieldinfo['shr06']  =  { 'fname'  : ['uzonal_6km','umeridional_6km','uzonal_surface','umeridional_surface'], 'filename': 'diag', 'skip':50 }
fieldinfo['shr01']  =  { 'fname'  : ['uzonal_1km','umeridional_1km','uzonal_surface','umeridional_surface'], 'filename': 'diag', 'skip':50 }

# Enter wind barb info for list of pressure levels
for plev in ['200', '250', '300', '500', '700', '850', '925']:
    fieldinfo['wind'+plev] = { 'fname' : ['uzonal_'+plev+'hPa', 'umeridional_'+plev+'hPa'], 'filename':'diag', 'skip':50}


def makeEnsembleListMPAS(wrfinit, timerange, ENS_SIZE, g193=False):
    # create lists of files (and missing file indices) for various file types
    shr, ehr = timerange
    file_list    = { 'wrfout':[], 'diag':[] }
    missing_list = { 'wrfout':[], 'diag':[] }
    missing_index = 0
    for hr in range(shr,ehr+1):
        wrfvalidstr = (wrfinit + datetime.timedelta(hours=hr)).strftime('%Y-%m-%d_%H.%M.%S')
        yyyymmddhh = wrfinit.strftime('%Y%m%d%H')
        for mem in range(1,ENS_SIZE+1):
            diag   = '/glade/p/nsc/nmmm0046/schwartz/MPAS_ens_15-3km_mesh/POST/%s/ens_%d/diag.%s.nc'%(yyyymmddhh,mem,wrfvalidstr)
            if g193:
                diag   = '/glade/p/nsc/nmmm0046/schwartz/MPAS_ens_15-3km_mesh/POST/%s/ens_%d/diag_latlon_g193.%s.nc'%(yyyymmddhh,mem,wrfvalidstr)
            logging.debug(diag)
            if os.path.exists(diag): file_list['diag'].append(diag)
            else: 
                missing_list['diag'].append(missing_index)
                logging.warning(f"{diag} does not exist")
            missing_index += 1
    logging.debug(f"file_list {file_list}")
    if not file_list['diag']:
        logging.info('Empty file_list')
    return (file_list, missing_list)
