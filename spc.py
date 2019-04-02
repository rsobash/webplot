import pandas as pd
import pytz
import glob
import sys # for stderr output
import pdb
import numpy as np
import datetime
import os # for basename
import matplotlib.path as mpath
import cartopy
import sqlite3

def spc_event_filename(event_type):
    name = '/glade/work/ahijevyc/share/SPC/'+event_type+'/'
    if event_type=='torn':
        return "/glade/work/ahijevyc/share/SPC/1950-2017_actual_tornadoes.csv"
    name = name + "1955-2017_"+event_type+".csv"
    return name



def get_storm_reports( 
        start = datetime.datetime(2016,6,10,0,0,0,0,pytz.UTC), 
        end   = datetime.datetime(2016,7,1,0,0,0,0,pytz.UTC), 
        event_types = ["torn", "wind", "hail"],
        debug      = False
        ):

    if debug:
        print("storm_reports: start:",start)
        print("storm_reports: end:",end)
        print("storm_reports: event types:",event_types)

    # Create one DataFrame to hold all event types.
    all_rpts = pd.DataFrame()
    for event_type in event_types:
        rpts_file = spc_event_filename(event_type)

        # csv format described in http://www.spc.noaa.gov/wcm/data/SPC_severe_database_description.pdf
        # SPC storm report files downloaded from http://www.spc.noaa.gov/wcm/#data to 
        # cheyenne:/glade/work/ahijevyc/share/ Mar 2019.
        # Multi-year zip files have headers; pre-2016 single-year csv files have no header. 

        dtype = {
                "om": np.int64,
                "yr": np.int32,
                "mo": np.int32,
                "dy": np.int32,
                "date": str,
                "tz": np.int32,
                "st": str,
                "stf": np.int64, # State FIPS number. some Puerto Rico codes are incorrect
                "stn": np.int64,  # State number - number of this tornado in this state in this year
                "mag": np.float64, # you might think there is "sz" for hail and "f" for torn, but all "mag"
                "inj": np.int64,
                "fat": np.int64,
                "loss": np.float64,
                "closs": np.float64,
                "slat": np.float64,
                "slon": np.float64,
                "elat": np.float64,
                "elon": np.float64,
                "len": np.float64,
                "wid": np.float64,
                "ns": np.int64,
                "sn": np.int64,
                "sg": np.int64,
                "f1": np.int64,
                "f2": np.int64,
                "f3": np.int64,
                "f4": np.int64,
                "mt": str,
                "fc": np.int64,
                }
        rpts = pd.read_csv(rpts_file, parse_dates=[['date','time']], dtype=dtype, infer_datetime_format=True)
        if debug:
            print("input file:",rpts_file)
            print("read",len(rpts),"lines")
            pdb.set_trace()
        rpts["event_type"] = event_type
        rpts["source"] = os.path.basename(rpts_file)

        # -9 = unknown tornado F-scale
        # Change -9 to NaN
        rpts["mag"].replace(to_replace=-9, value = np.nan, inplace=True)

        # Augment event type for large hail and high wind.
        if event_type == "hail":
            largehail = rpts.mag >= 2.
            if any(largehail):
                rpts.loc[largehail,"event_type"] = "large hail"
        if event_type == "wind":
            highwind = rpts.mag >= 65.
            if any(highwind):
                rpts.loc[highwind, "event_type"] = "high wind"

        # Fix tz=6.  This is MDT in the NECI Storm Events database.
        # This is confusing. A time in MDT is the same as the time in CST.
        # So simply change tz to 3 (CST).
        MDT = rpts['tz'] == 6
        if any(MDT):
            MDT_timezones = rpts[['om','date_time','tz','event_type', 'source']][MDT]
            print("get_storm_reports: found",len(MDT_timezones),"MDT time zones")
            if debug:
                print(MDT_timezones, file=sys.stderr)
                print("changing tz from 6 to 3 because CST=MDT")
                print("WARNING - proceeding with program. Wrote to SPC Apr 1 2019 about fixing these lines", file=sys.stderr)
            Email20190401PatrickMarsh = "...took over database in 2017 and have no record or documentation as to what those timezones are. Each year I append new information to the end of the old information, so timezones will continue to exist as is until I can learn what those time zones are.  take a look at the NCEI version of storm data. They may have information I do not."
            rpts.loc[MDT == 6, "tz"] = 3
        # When timezone = 0, it is 'UNK' (unknown?) in the NCEI Storm Events database.
        # When timezone = 6, it is 'MDT' (Mountain Daylight Time?) in the NCEI Storm Events database.
        """
        om       yr_mo_dy_time  tz county/zone time zone according to NCEI Storm Events
2507   245 1956-06-01 11:33:00   0 UNK
2855    89 1957-04-02 23:45:00   0 UNK
8145   216 1965-05-05 14:45:00   6 MDT
9569   158 1967-04-21 12:33:00   0 UNK
13409  264 1972-05-13 18:08:00   0 UNK
15792  804 1974-08-13 15:03:00   0 UNK
20640  458 1980-06-04 16:30:00   0 UNK
22101  271 1982-05-11 14:25:00   0 UNK
24007  200 1984-04-26 19:32:00   0 CST
25899  501 1986-07-01 22:15:00   6 MDT
26054  656 1986-09-04 18:55:00   6 MDT
28793  419 1990-05-24 15:00:00   6 MDT
28797  422 1990-05-24 16:00:00   6 MDT
28810  433 1990-05-24 18:33:00   6 MDT
29192  815 1990-06-27 20:00:00   6 MDT
29232  855 1990-07-05 21:10:00   6 MDT
29347  971 1990-08-15 18:30:00   6 MDT
33262  151 1994-04-22 18:06:00   6 MDT
33557  446 1994-05-31 15:00:00   6 MDT
33572  461 1994-06-06 14:40:00   6 MDT
33574  463 1994-06-06 15:00:00   6 MDT
33585  474 1994-06-07 14:47:00   6 MDT
33586  476 1994-06-07 15:57:00   6 MDT
33587  477 1994-06-07 16:10:00   6 MDT
33589  478 1994-06-07 16:35:00   6 MDT
33592  482 1994-06-07 18:50:00   6 not in Storm Events database
33783  672 1994-06-29 15:45:00   6 MDT
33834  722 1994-07-06 18:17:00   6 MDT
33855  744 1994-07-12 14:30:00   6 MDT
33886  775 1994-07-18 15:30:00   6 MDT
33887  776 1994-07-18 16:00:00   6 MDT
33888  777 1994-07-18 16:25:00   6 MDT
33889  778 1994-07-18 16:40:00   6 MDT
33890  779 1994-07-18 16:55:00   6 not in Storm Events database
33891  780 1994-07-18 16:55:00   6 not in Storm Events database
33892  781 1994-07-18 17:00:00   6 MDT
"""            


        # All times, except for ?=unknown and 9=GMT, were converted to 3=CST.
        # tz=3 is CST
        # tz=9 is GMT
        # Convert to UTC by adding 9 and subtracting tz hours.
        rpts["time"] = rpts['date_time'] + pd.to_timedelta(9 - rpts.tz, unit='h')
        # make time-zone aware datetime object
        rpts["time"] = rpts["time"].dt.tz_localize(pytz.UTC)

        if any(rpts['tz'] == 0):
            print("reports file: "+rpts_file, file=sys.stderr)
            unknown_timezones = rpts[['om','date_time','tz','event_type', 'source']][rpts['tz'] == 0]
            if debug:
                print("get_storm_reports: found",len(unknown_timezones),"unknown time zones")
                print(unknown_timezones, file=sys.stderr)

        time_window = (rpts.time >= start) & (rpts.time < end)
        rpts = rpts[time_window]
        if debug:
            print("found",len(rpts),event_type,"reports")

        # Sanity Check:
        # Verify I get the same thing as Ryan Sobash's sqlite3 database
        conn = sqlite3.connect("/glade/u/home/sobash/2013RT/REPORTS/reports_all.db")
        sqltable = "reports_" + event_type
        # Could apply a datetime range here (WHERE datetime BETWEEN yyyy/mm/dd hh:mm:ss and blah), but converting from UTC to CST is tricky.
        sqlcommand = "SELECT * FROM "+sqltable
        sql_df = pd.read_sql_query(sqlcommand, conn, parse_dates=['datetime'])
        conn.close()
        # Add 6 hours to datetime. This converts to UTC.
        sql_df["datetime"] = sql_df["datetime"] + pd.to_timedelta(6, unit='h')
        # make it aware of its UTC timezone.
        sql_df["datetime"] = sql_df["datetime"].dt.tz_localize(pytz.UTC)
        sql_df = sql_df[(sql_df.datetime >= start) & (sql_df.datetime < end)]

        # See if they have the same number of rows
        if len(sql_df) != len(rpts):
            print("My data don't match Ryan's SQL database")
            pdb.set_trace()
            sys.exit(1)
        # See if they have the same times
        if (sql_df["datetime"].values != rpts["time"].values).any():
            print("My data times don't match Ryan's SQL database")
        # See if they have the same locations
        same_columns = ["slat", "slon", "elat", "elon"]
        if (sql_df[same_columns].values != rpts[same_columns].values).any():
            print("My data locations don't match Ryan's SQL database")
            max_abs_difference = np.max(np.abs(sql_df[same_columns]-rpts[same_columns]))
            print("max abs difference")
            print(max_abs_difference)
            if all(max_abs_difference < 0.000001):
                print("who cares about such a small difference?")
            else:
                pdb.set_trace()
                sys.exit(1)
        if debug:
            print("From Ryan's SQL database")
            print(sqlcommand)
            print(sql_df.to_string())


        all_rpts = all_rpts.append(rpts, ignore_index=True, sort=False)

    return all_rpts

def plot(storm_reports, ax, scale=1, drawradius=0, alpha=0.4, debug=False):

    if storm_reports.empty:
        # is this the right thing to return? what about empty list []? or rpts?
        if debug:
            print("spc.plot(): storm reports DataFrame is empty. Returning")
        return

    # Color, size, marker, and label of wind, hail, and tornado storm reports
    kwdict = {
            "wind":       {"c" : 'blue',  "s": 8*scale, "marker":"s", "label":"Wind Report"},
            "high wind":  {"c" : 'black', "s":12*scale, "marker":"s", "label":"Wind Report/HI"},
            "hail":       {"c" : 'green', "s":12*scale, "marker":"^", "label":"Hail Report"},
            "large hail": {"c" : 'black', "s":16*scale, "marker":"^", "label":"Hail Report/LG"},
            "torn":       {"c" : 'red',   "s":12*scale, "marker":"v", "label":"Torn Report"}
            }

    storm_rpts_plots = []

    for event_type in ["wind", "high wind", "hail", "large hail", "torn"]:
        if debug:
            print("looking for "+event_type)
        xrpts = storm_reports[storm_reports.event_type == event_type]
        print("plot",len(xrpts),event_type,"reports")
        if len(xrpts) == 0:
            continue
        lons, lats = xrpts.slon.values, xrpts.slat.values
        storm_rpts_plot = ax.scatter(lons, lats, alpha = alpha, edgecolors="None", **kwdict[event_type],
                transform=cartopy.crs.Geodetic())
        storm_rpts_plots.append(storm_rpts_plot)
        if debug:
            pdb.set_trace()
        if drawradius > 0:
            # With lons and lats, specifying more than one dimension allows individual points to be drawn. 
            # Otherwise a grid of circles will be drawn.
            # It warns about using PlateCarree to approximate Geodetic. It still warps the circles
            # appropriately, so I think this is okay.
            within_radius = ax.tissot(rad_km = drawradius, lons=lons[np.newaxis], lats=lats[np.newaxis], 
                facecolor=kwdict[event_type]["c"], alpha=0.4, label=str(drawradius)+" km radius")
            # TODO: Legend does not support tissot cartopy.mpl.feature_artist.
            # A proxy artist may be used instead.
            # matplotlib.org/users/legend_guide.html#
            # creating-artists-specifically-for-adding-to-the-legend-aka-proxy-artists
            # storm_rpts_plots.append(within_radius)

    return storm_rpts_plots



def to_MET(df, gribcode=187):
    # gribcode 187 :lightning
    # INPUT: storm_reports DataFrame from spc.get_storm_reports()
    # Output: MET point observation format, one observation per line.
    # Use gribcode if specified. Otherwise default.
    # 10 = WIND
    # 
    # Each observation line will consist of the following 11 columns of data:
    met_columns = ["Message_Type", "Station_ID", "Valid_Time", "Lat", "Lon", "Elevation", "Grib_Code", "Level", "Height", "QC_String", "Observation_Value"]
    
    df["Message_Type"] = "ADPSFC"
    df["Station_ID"]  = df.st
    df["Valid_Time"] = df.time.dt.strftime('%Y%m%d_%H%M%S') #df.time should be aware of UTC time zone
    df["Lat"]  = (df.slat+df.elat)/2
    df["Lon"]  = (df.slon+df.elon)/2
    df["Elevation"]  = "NA"
    df["Grib_Code"] = gribcode
    df["Height"] = "NA"
    # Change grib_code to 10 (WIND) where event type is wind
    df.loc[df['event_type'].str.contains('wind'), "Grib_Code"] = 10
    df.loc[df['event_type'].str.contains('wind'), "Height"] = 10
    df["Level"] = "NA"
    df["QC_String"] = "NA"
    df["Observation_Value"] = df.mag
    # index=False don't write index number
    # Change NaN to "NA" MET considers "NA" missing
    # This may not matter, but by adding a string 'NA', it changes the format of the entire column. Floats change to integers (or maybe strings). 
    df.replace(to_replace=np.nan, value = 'NA', inplace=True)
    return df.to_string(columns=met_columns,index=False,header=False) + "\n"



# NCDC storm event details don't have a lat and lon. Just a city, a range, and an azimuth.
# Perhaps the lat and lon are in the storm event location files (as opposed to "details").
# Maybe just go with SPC storm reports.

def stormEvents(
        idir="/glade/work/ahijevyc/share/stormevents", 
        start = datetime.datetime(2016,6,1), 
        end   = datetime.datetime(2016,7,1), 
        event_type = "Tornado", 
        debug      = False
        ):
    # INPUT: idir
    # Directory with input files. For example: StormEvents_details-ftp_v1.0_d2015_c20180525.csv
    # Downloaded from NCDC storm events database at https://www.ncdc.noaa.gov/stormevents/ftp.jsp

    ifiles = glob.glob(idir+"/StormEvents_details-*csv*")
    year = start.year
    ifile = [s for s in ifiles if "_d"+str(year)+"_c" in s]
    if len(ifile) != 1:
        print("Did not find one storm events file", ifile)
        pdb.set_trace()

    ifile = ifile[0]
    if debug:
        pdb.set_trace()

    events = pd.read_csv(ifile)
    # This is ugly
    events["BEGIN"] = pd.to_datetime(events.BEGIN_DATE_TIME, format="%d-%b-%y %H:%M:%S")
    events["END"]  = pd.to_datetime(events.END_DATE_TIME, format="%d-%b-%y %H:%M:%S")
    # Just return a particular event type. 'Thunderstorm Wind' or 'Tornado' or 'Hail', for example.
    if event_type is not None:
        events = events[events.EVENT_TYPE == event_type]

    if start is not None:
        events = events[events.BEGIN >= start]

    if end is not None:
        events = events[events.END < end]

    return events

def events2met(df):
    # Point Observation Format
    # input ASCII MET point observation format contains one observation per line. 
    # Each input observation line should consist of the following 11 columns of data:
    met_columns = ["Message_Type", "Station_ID", "Valid_Time", "Lat", "Lon", "Elevation", "Grib_Code", "Level", "Height", "QC_String", "Observation_Value"]
    
    df["Message_Type"] = "ADPSFC"
    df["Station_ID"]  = df.WFO
    df["Valid_Time"] = df.BEGIN.dt.strftime('%Y%m%d_%H%M%S')
    df["Lat"]  = df.lat # I don't know how to get this. 
    df["Lon"]  = df.lon
    df["Elevation"]  = "NA"
    df["Grib_Code"] = "NA"
    df["Level"] = "NA"
    df["Height"] = "NA"
    df["QC_String"] = "NA"
    df["Observation_Value"] = df.TOR_F_SCALE
    print(df.to_string(columns=met_columns))


