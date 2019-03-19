import pandas as pd
import pdb
import re
import csv
import os, sys
from netCDF4 import Dataset
import numpy as np


def kts2category(kts):
   category = -1
   if kts > 34:
      category = 0
   if kts > 64:
      category = 1
   if kts > 83:
      category = 2
   if kts > 96:
      category = 3
   if kts > 113:
      category = 4
   if kts > 137:
      category = 5
   return category


ifile = '/glade/work/ahijevyc/work/atcf/Irma.ECMWF.dat'
ifile = '/glade/scratch/mpasrt/uni/2018071700/latlon_0.500deg_0.25km/gfdl_tracker/tcgen/fort.64'



# Standard ATCF columns (doesn't include track id, like in fort.66).
# https://www.nrlmry.navy.mil/atcf_web/docs/database/new/abrdeck.html
atcfcolumns=["basin","cy","initial_time","technum","model","fhr","lat","lon","vmax","minp","ty",
    "rad", "windcode", "rad1", "rad2", "rad3", "rad4", "pouter", "router", "rmw", "gusts", "eye",
    "subregion", "maxseas", "initials", "dir", "speed", "stormname", "depth", "seas", "seascode",
    "seas1", "seas2", "seas3", "seas4", "userdefined", "userdata"]

def cyclone_phase_space_columns():
    names = []
    names.append('cpsB') # Cyclone Phase Space "Parameter B" for thermal asymmetry. (Values are *10)
    names.append('cpsll') # Cyclone Phase Space lower level (600-900 mb) thermal wind parameter, for diagnosing low-level warm core. (Values are *10)
    names.append('cpsul') # Cyclone Phase Space upper level (300-600 mb) thermal wind parameter, for diagnosing upper-level warm core. (Values are *10)
    return names 


def getcy(cys):
    return cys[0:2]

def read(ifile = ifile, debug=False, fullcircle=False):
    # Read data into Pandas Dataframe
    print('reading', ifile, 'fullcircle=', fullcircle)
    names = list(atcfcolumns) # make a copy of list, not a copy of the reference to the list.
    converters={
            # The problem with CY is ATCF only reserves 2 characters for it.
            "cy" : lambda x: x.strip(), # cy is not always an integer (e.g. 10E) # Why strip leading zeros?
            "initial_time" : lambda x: pd.to_datetime(x.strip(),format='%Y%m%d%H'),
            "technum" : lambda x: x.strip(), #strip leading spaces but not leading zeros
            "model" : lambda x: x.strip(), # Strip leading whitespace - for matching later.
            "vmax": float,
            "minp": float,
            "minp": float,
            "ty" : lambda x: x.strip(),
            "windcode" : lambda x: x[-3:],
            "rad1": float,
            "rad2": float,
            "rad3": float,
            "rad4": float,
            "pouter": float,
            "router": float,
            "subregion": lambda x: x[-2:],
             # subregion ends up being 3 when written with .to_string
             # strange subregion only needs one character, but official a-decks leave 3. 
            "initials" : lambda x: x[-3:],
            'stormname': lambda x: x[-9:],
            'depth'    : lambda x: x[-1:],
            "seascode" : lambda x: x[-3:],
            "seas1": float,
            "seas2": float,
            "seas3": float,
            "seas4": float,
            "userdefined" : lambda x: x.strip(),
            "userdata" : lambda x: x.strip(),
        }
    dtype={
            'rmw'      : np.float64,
            'gusts'    : np.float64,
            'eye'      : np.float64,
            'maxseas'  : np.float64,
            'dir'      : np.float64,
            'speed'    : np.float64,
            "seas"     : np.float64,
           } 

    # Tried using converter for these columns, but couldn't convert 4-space string to float.
    # If you add a key na_values, also add it to dtype dict, and remove it from converters.
    na_values = {
            "rmw"     : 4*' ',
            "gusts"   : 4*' ',
            "eye"     : 4*' ',
            "maxseas" : 4*' ',
            "dir"     : 4*' ',
            "speed"   : 4*' ',
            "seas"    : 3*' ', # one less than other columns
            }


    reader = csv.reader(open(ifile),delimiter=',')
    testline = next(reader)
    num_cols = len(testline)
    if debug:
        print(testline)
        print('num_cols=',num_cols)
    del reader

    # Output from HWRF vortex tracker, fort.64 and fort.66
    # are mostly ATCF format but have subset of columns
    if num_cols == 43:
        print('assume HWRF tracker fort.64-style output with 43 columns in', ifile)
        TPstr = "THERMO PARAMS"
        if testline[35].strip() != TPstr:
            print("expected 36th column to be", TPstr)
            print("got", testline[35].strip())
            sys.exit(4)
        for ii in range(20,35):
            names[ii] = "space filler"
        names = names[0:35]
        names.append(TPstr)
        names.extend(cyclone_phase_space_columns())
        names.append('warmcore')
        names.append("warmcore_strength")
        names.append("string")
        names.append("string")

    # fort.66 has track id in the 3rd column.
    if num_cols == 31:
        print('assume fort.66-style with 31 columns in', ifile)
        # There is a cyclogenesis ID column for fort.66
        if debug:
            print('inserted ID for cyclogenesis in column 2')
        names.insert(2, 'id') # ID for the cyclogenesis
        print('Using 1st 21 elements of names list')
        names = names[0:21]
        if debug:
            print('redefining columns 22-31')
        names.extend(cyclone_phase_space_columns())
        names.append('warmcore')
        names.append('dir')
        names.append('speedms')
        names.append('vort850mb')
        names.append('maxvort850mb')
        names.append('vort700mb')
        names.append('maxvort700mb')

    # TODO read IDL output
    if num_cols == 44 and 'min_warmcore_fract d' in testline[35]:
        if debug:
            print("Looks like IDL output")
        names = [n.replace('userdata', 'min_warmcore_fract') for n in names]
        names.append('dT500')
        names.append('dT200')
        names.append('ddZ850200')
        names.append('rainc')
        names.append('rainnc')
        names.append('id')

    if num_cols == 11:
        print("Assuming simple adeck with 11 columns")
        if ifile[-4:] != '.dat':
            print("even though file doesn't end in .dat", ifile)
        names = names[0:11]

    usecols = list(range(len(names)))

    # If you get a beyond index range (or something like that) error, see if userdata column is intermittent and has commas in it. 
    # If so, clean it up (i.e. truncate it)

    #atrack = ['basin','cy','initial_time','technum','model']
    #if 'id' in names:
    #    atrack.append('id')

    if debug:
        print("before pd.read_csv")
        print('column names', names)
        print("usecols=",usecols)
        print("converters=",converters)
        print("dype=", dtype)
    df = pd.read_csv(ifile,index_col=None,header=None, delimiter=",", usecols=usecols, names=names, 
            converters=converters, na_values=na_values, dtype=dtype) 
    # fort.64 has asterisks sometimes. Problem with hwrf_tracker. 
    badlines = df['lon'].str.contains("\*")
    if any(badlines):
        df = df[~badlines]

    # Extract last character of lat and lon columns
    # Multiply integer by -1 if "S" or "W"
    # Divide by 10
    S = df.lat.str[-1] == 'S'
    lat = df.lat.str[:-1].astype(float) / 10.
    lat[S] = lat[S] * -1
    df.lat = lat
    W = df.lon.str[-1] == 'W'
    lon = df.lon.str[:-1].astype(float) / 10.
    lon[W] = lon[W] * -1
    df.lon = lon

    if debug:
        pdb.set_trace()
    # Derive valid time.   valid_time = initial_time + fhr
    # Use datetime module to add, where yyyymmddh is a datetime object and fhr is a timedelta object.
    df['valid_time'] = df.initial_time + pd.to_timedelta(df.fhr, unit='h')

    for col in atcfcolumns:
        if col not in df.columns:
            if debug:
                print(col, 'not in DataFrame. Fill with appropriate value.')
            # if 'rad' column doesn't exist make it zeroes
            if col in ['rad', 'rad1', 'rad2', 'rad3', 'rad4','pouter', 'router', 'seas', 'seas1','seas2','seas3','seas4']:
                df[col] = 0.

            # Initialize other default values.
            if col in ['windcode', 'seascode']:
                df[col] = '   '

            # Numbers are NaN
            if col in ['rmw','gusts','eye','maxseas','dir','speed']:
                df[col] = np.NaN

            # Strings are empty
            if col in ['subregion','stormname','depth','userdefined','userdata']:
                df[col] = ''

            if col in ['initials', 'depth']:
                df[col] = 'X'
            

    if debug:
        pdb.set_trace()
    if fullcircle:
        if debug:
            print("full circle wind radii")
        # Full circle wind radii instead of quadrants
        # TODO better way than this hack
        df['windcode'] = 'AAA'
        df['rad1'] = df[['rad1','rad2','rad3','rad4']].max(axis=1)
        df['rad2'] = 0
        df['rad3'] = 0
        df['rad4'] = 0

    return df

def x2s(x):
    # Convert absolute value of float to integer number of tenths for ATCF lat/lon
    # called by lat2s and lon2s
    x *= 10
    x = np.around(x)
    x = np.abs(x)
    return str(int(x))

def lat2s(lat):
    NS = 'N' if lat >= 0 else 'S'
    lat = x2s(lat) + NS + ','
    return lat

def lon2s(lon):
    EW = 'E' if lon >= 0 else 'W'
    lon = '%4s' % x2s(lon) + EW + ','
    return lon

# function to compute great circle distance between point lat1 and lon1 and arrays of points 
# given by lons, lats
# Returns 2 things:
#   1) distance in km
#   2) initial bearing from 1st pt to 2nd pt.
def dist_bearing(lon1,lons,lat1,lats):
    lon1 = np.radians(lon1)
    lons = np.radians(lons)
    lat1 = np.radians(lat1)
    lats = np.radians(lats)
    # great circle distance. 
    arg = np.sin(lat1)*np.sin(lats)+np.cos(lat1)*np.cos(lats)*np.cos(lon1-lons) 
    #arg = np.where(np.fabs(arg) < 1., arg, 0.999999) 

    dlon = lons-lon1
    bearing = np.arctan2(np.sin(dlon)*np.cos(lats), np.cos(lat1)*np.sin(lats) - np.sin(lat1)*np.cos(lats)*np.cos(dlon)) 

    # convert from radians to degrees
    bearing = np.degrees(bearing)

    # -180 - 180 -> 0 - 360
    bearing = (bearing + 360) % 360
    
    # Ellipsoid [CLARKE 1866]  Semi-Major Axis (Equatorial Radius)
    a = 6378.2064
    return np.arccos(arg)* a, bearing 



ms2kts = 1.94384
km2nm = 0.539957

quads = {'NE':0, 'SE':90, 'SW':180, 'NW':270}
thresh_kts = np.array([34, 50, 64])



def get_max_ext_of_wind(speed_kts, distance_km, bearing, raw_vmax_kts, quads=quads, thresh_kts=thresh_kts, debug=False):
    rad_nm = {}
    # Put in dictionary "rad_nm" where rad_nm = {
    #                                       34: {'NE':rad1, 'SE':rad2, 'SW':rad3, 'NW':rad4},
    #                                       50: {'NE':rad1, 'SE':rad2, 'SW':rad3, 'NW':rad4},
    #                                       64: {'NE':rad1, 'SE':rad2, 'SW':rad3, 'NW':rad4} 
    #                                     }
    rad_nm['raw_vmax_kts'] = raw_vmax_kts
    rad_nm['thresh_kts'] = thresh_kts
    rad_nm['quads'] = quads
    for wind_thresh_kts in thresh_kts[thresh_kts < raw_vmax_kts]:
        rad_nm[wind_thresh_kts] = {}
        for quad,az in quads.items():
            # Originally had distance_km < 800, but Chris D. suggested 300nm in Sep 2018 email
            # This was to deal with Irma and the unrelated 34 knot onshore flow in Georgia
            # Looking at HURDAT2 R34 sizes (since 2004), ex-tropical storm Karen 2015 had 710nm.
            # Removing EX storms, the max was 480 nm in Hurricane Sandy 2012
            iquad = (az <= bearing) & (bearing < az+90) & (speed_kts >= wind_thresh_kts) & (distance_km < 300./km2nm)
            rad_nm[wind_thresh_kts][quad] = 0
            if np.sum(iquad) > 0:
                max_dist_of_wind_threshold = np.max(distance_km[iquad]) * km2nm
                imax_dist_of_wind_threshold = np.argmax(distance_km[iquad])
                if debug:
                    print('get_max_ext_of_wind():', wind_thresh_kts, quad, '%3d-%3d'%(az,az+90), '%4d'%np.sum(iquad), '%10.6f'%max_dist_of_wind_threshold, '%10.6f'%bearing[iquad][imax_dist_of_wind_threshold])
                rad_nm[wind_thresh_kts][quad] = max_dist_of_wind_threshold
    return rad_nm


def derived_winds(u10, v10, mslp, lonCell, latCell, row, vmax_search_radius=250., mslp_search_radius=100., debug=False):


    # Given a row (with row.lon and row.lat)...

    # Derive cell distances and bearings
    distance_km, bearing = dist_bearing(row.lon, lonCell, row.lat, latCell)

    # Derive 10m wind speed and Vt from u10 and v10
    speed_kts = np.sqrt(u10**2 + v10**2) * ms2kts

    # Tangential (cyclonic) wind speed
    # v dx - u dy 
    dx = lonCell - row.lon
    # work on the dateline?
    dx[dx>=180] = dx[dx>=180]-360.
    dy = latCell - row.lat
    Vt = v10 * dx - u10 * dy
    if row.lat < 0:
        Vt = -Vt

    # Restrict Vmax search
    vmaxrad = distance_km < vmax_search_radius
    ispeed_max = np.argmax(speed_kts[vmaxrad])
    raw_vmax_kts =  speed_kts[vmaxrad].max()

    # Check if tangential component of max wind is negative (anti-cyclonic)
    if Vt[vmaxrad][ispeed_max] < 0:
        print("center", row.valid_time, row.lat, row.lon)
        print("max wind is anti-cyclonic!", Vt[vmaxrad][ispeed_max])
        print("max wind lat/lon", latCell[vmaxrad][ispeed_max], lonCell[vmaxrad][ispeed_max])
        print("max wind U/V",         u10[vmaxrad][ispeed_max],     v10[vmaxrad][ispeed_max])
        if debug: pdb.set_trace()

    # Check if average tangential wind within search radius is negative (anti-cyclonic)
    average_tangential_wind = np.average(Vt[vmaxrad])
    if average_tangential_wind < 0:
        print("center", row.valid_time, row.lat, row.lon)
        print("avg wind is anti-cyclonic!", average_tangential_wind)
        if debug: pdb.set_trace()

    # Get radius of max wind
    raw_RMW_nm = distance_km[vmaxrad][ispeed_max] * km2nm
    if debug:
        print('max wind lat', latCell[vmaxrad][ispeed_max], 'lon', lonCell[vmaxrad][ispeed_max])

    # Restrict min mslp search
    mslprad = distance_km < mslp_search_radius
    raw_minp = mslp[mslprad].min() / 100.

    # Get max extent of wind at thresh_kts thresholds.
    rad_nm = get_max_ext_of_wind(speed_kts, distance_km, bearing, raw_vmax_kts, debug=debug)

    return raw_vmax_kts, raw_RMW_nm, raw_minp, rad_nm


def add_wind_rad_lines(row, rad_nm, fullcircle=False, debug=False):
    raw_vmax_kts = rad_nm['raw_vmax_kts']
    thresh_kts = rad_nm['thresh_kts']
    # if not empty...must be NEQ
    if row.windcode.strip() and row.windcode != 'NEQ':
        print('bad windcode', row.windcode, 'in', row)
        print('expected NEQ')
        sys.exit(1)
    lines = pd.DataFrame()
    for thresh in thresh_kts[thresh_kts < raw_vmax_kts]:
        if any(rad_nm[thresh].values()):
            newrow = row.copy()
            newrow[['windcode','rad','rad1','rad2','rad3','rad4']] = ['NEQ', thresh, rad_nm[thresh]['NE'], rad_nm[thresh]['SE'], rad_nm[thresh]['SW'], rad_nm[thresh]['NW']]
            # Append row with 34, 50, or 64 knot radii
            lines = lines.append(newrow)
            if fullcircle:
                # Append row with full circle 34, 50, or 64 knot radius
                # MET-TC will not derive this on its own - see email from John Halley-Gotway Oct 11, 2018
                # Probably shouldn't have AAA and NEQ in same file. 
                newrow = row.copy()
                newrow[['windcode','rad','rad1']] = ['AAA',thresh, np.nanmax(list(rad_nm[thresh].values()))] 
                lines = lines.append(newrow)
    
    return lines


def update_df(df, row, raw_vmax_kts, raw_RMW_nm, raw_minp, rad_nm, debug=False):

    # Called by origgrid and origmesh

    if debug:
        print('before', row[['valid_time','lon','lat', 'vmax', 'minp', 'rmw']]) 
    row.vmax = raw_vmax_kts
    row.minp = raw_minp
    row.rmw  = raw_RMW_nm
    # Add note of original mesh = True in user data (not defined) column
    row.userdata += 'origmeshTrue'
    if debug:
        print('after', row[['vmax', 'minp', 'rmw']]) 

    # Can't get rid of SettingWithCopyWarning! 
    df.loc[row.name,:] = row


    # Append 34/50/64 knot lines to DataFrame
    newlines = add_wind_rad_lines(row, rad_nm, debug=debug)
    # If there are new lines, drop the old one and append new ones. 
    if not newlines.empty:
        df.drop(row.name, inplace=True)
        df = df.append(newlines)



    # Sort DataFrame by index (deal with appended wind radii lines)
    # sort by rad too
    df = df.sort_index().sort_values(['initial_time','fhr','rad'])

    return df





def write(ofile, df, fullcircle=False, debug=False):
    if df.empty:
        print("afcf.write(): DataFrame is empty.", ofile, "not written")
        return

    # TODO: deal with fullcircle.
    print("writing", ofile)

    # Valid time is not part of ATCF file.
    del(df["valid_time"])

    if debug:
        pdb.set_trace()
    # Had to add parentheses () to .len to not get error about instancemethod not being iterable.
    if max(df.cy.str.len()) > 2:
        print('cy more than 2 characters. Truncating...')
        # Append full CY to userdata column (after a comma)
        df['userdata'] = df['userdata'] + ', ' + df['cy']
        # Keep first 2 characters
        df['cy'] = df['cy'].str.slice(0,2)
    formatters={
            "basin": '{},'.format,
            # The problem with CY is ATCF only reserves 2 characters for it.
            "cy": lambda x: x.zfill(2)+"," , # not always an integer (e.g. 10E) # 20181116 force to be integer
            # Convert initial_time from datetime to string.
            "initial_time":lambda x: x.strftime('%Y%m%d%H,'),
            "technum":'{},'.format,
            "model":'{},'.format,
            "fhr":'{:3.0f},'.format,
            "lat":lat2s, 
            "lon":lon2s, 
            "vmax":'{:3.0f},'.format,
            "minp":'{:4.0f},'.format,
            "ty":'{},'.format,
            "windcode":'{:>3s},'.format,
            "rad":'{:3.0f},'.format,
            "rad1":'{:4.0f},'.format,
            "rad2":'{:4.0f},'.format,
            "rad3":'{:4.0f},'.format,
            "rad4":'{:4.0f},'.format,
            "pouter":'{:4.0f},'.format,
            "router":'{:4.0f},'.format,
            "rmw":'{:3.0f},'.format,
            "gusts":'{:3.0f},'.format,
            "eye":'{:3.0f},'.format,
            "subregion":'{:>2s},'.format,
            "maxseas":'{:3.0f},'.format,
            "initials":'{:>3s},'.format,
            "dir":'{:3.0f},'.format,
            "speed":'{:3.0f},'.format,
            "stormname":'{:>9s},'.format,
            "depth":'{:>1s},'.format,
            "seas":'{:2.0f},'.format,
            "seascode":'{:>3s},'.format,
            "seas1":'{:4.0f},'.format,
            "seas2":'{:4.0f},'.format,
            "seas3":'{:4.0f},'.format,
            "seas4":'{:4.0f},'.format,
            "userdefined":'{:>18s},'.format,
            "cpsB":'{:4.0f},'.format,
            "cpsll":'{:4.0f},'.format,
            "cpsul":'{:4.0f},'.format,
            "warmcore":'{},'.format,
            "direction":'{},'.format,
            }
    junk = df.to_string(header=False, index=False, na_rep=' ', columns=atcfcolumns, formatters=formatters)

    # TODO: FIX IDL-STYLE COLUMNS. STORMNAME IS MISSING A LETTER
    
    # na_rep=' ' has no effect
    # strings have extra space in front of them
    junk = junk.split('\n')
    # replace first 3 occurrences of ',  ' with ', '.
    # replace ',  XX,' with ', XX,'
    # replace 'nan' with '   '
    junk = [j.replace(',  ', ', ', 3).replace(',  XX,',', XX,').replace('nan','   ') for j in junk]
    #delete space before windcode
    junk = [j[:68]+j[69:] for j in junk]
    #delete space before initials
    junk = [j[:133]+j[134:] for j in junk]
    #delete space before depth
    junk = [j[:161]+j[162:] for j in junk]
    #delete space before seascode
    junk = [j[:168]+j[169:] for j in junk]
    junk = '\n'.join(junk)
    
    if debug:
        pdb.set_trace()

    f = open(ofile, "w")
    f.write(junk+"\n")

    f.close()
    print("wrote", ofile)


def origgrid(df, griddir, debug=False):
    # Get vmax, minp, radius of max wind, max radii of wind thresholds from ECMWF grid, not from tracker.
    # Assumes
    #   ECMWF data came from TIGGE and were converted from GRIB to netCDF with ncl_convert2nc.
    #   4-character model string in ATCF file is "EExx" (where xx is the 2-digit ensemble member).
    #   ECMWF ensemble member in directory named "ens_xx" (where xx is the 2-digit ensemble member). 
    #   File path is "ens_xx/${gs}yyyymmddhh.xx.nc", where ${gs} is the grid spacing (0p15, 0p25, or 0p5).

    for run_id, group in df.groupby(['initial_time', 'model']):
        initial_time, model = run_id
        m = re.search(r'EE(\d\d)', model)
        if not m:
            if debug:
                print('Assuming ECMWF ensemble member, but did not find EE\d\d in model string')
                print('no original grid for',model,'- skipping')
            continue
        ens = int(m.group(1)) # strip leading zero
        if ens < 1:
            continue
        # Allow some naming conventions
        # ens_n/yyyymmddhh.n.nc
        # ens_n/0p15yyyymmddhh_sfc.nc
        # ens_n/0p25yyyymmddhh_sfc.nc
        # ens_n/0p5yyyymmddhh_sfc.nc
        yyyymmddhh = initial_time.strftime('%Y%m%d%H')
        # If first filename doesn't exist, try the next one, and so on...
        # List in order of most preferred to least preferred.
        potential_gridfiles = [
                               "ens_"+str(ens)+"/"+ "0p15"+yyyymmddhh+"."+str(ens)+".nc",
                               "ens_"+str(ens)+"/"+ "0p25"+yyyymmddhh+"."+str(ens)+".nc",
                               "ens_"+str(ens)+"/"+ "0p5"+yyyymmddhh+"."+str(ens)+".nc",
                               "ens_"+str(ens)+"/"+ yyyymmddhh+"."+str(ens)+".nc"
                               ]
        for gridfile in potential_gridfiles:
            if os.path.isfile(griddir + gridfile):
                break
            else:
                print("no", griddir + gridfile)

        print('opening', gridfile)
        nc = Dataset(griddir + gridfile, "r")
        lon = nc.variables['lon_0'][:]
        lat = nc.variables['lat_0'][:]
        lonCell,latCell = np.meshgrid(lon, lat)
        u10s  = nc.variables['10u_P1_L103_GLL0'][:]
        v10s  = nc.variables['10v_P1_L103_GLL0'][:]
        mslps = nc.variables['msl_P1_L101_GLL0'][:]
        model_forecast_times = nc.variables['forecast_time0'][:]
        nc.close()
        for index, row in group.iterrows():
            if not any(model_forecast_times == row.fhr):
                print(row.fhr, 'not in model file')
                continue
            itime = np.argmax(model_forecast_times == row.fhr)
            u10  =  u10s[itime,:,:]
            v10  =  v10s[itime,:,:]
            mslp = mslps[itime,:,:]

            # Extract vmax, RMW, minp, and radii of wind thresholds
            raw_vmax_kts, raw_RMW_nm, raw_minp, rad_nm = derived_winds(u10, v10, mslp, lonCell, latCell, row, debug=debug)


            df = update_df(df, row, raw_vmax_kts, raw_RMW_nm, raw_minp, rad_nm, debug=debug)

    return df

if __name__ == "__main__":
    read()
