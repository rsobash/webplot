import pandas as pd
import pytz
import glob
import sys # for stderr output
import pdb
import numpy as np
import datetime
import matplotlib.path as mpath
import cartopy

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
        if debug:
            print("input file:",rpts_file)
            pdb.set_trace()

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
                "stf": np.int64,
                "stn": np.int64,
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
        rpts["event_type"] = event_type

        # Augment event type for large hail and high wind.
        if event_type == "hail":
            largehail = rpts.mag >= 2.
            if any(largehail):
                rpts.loc[largehail,"event_type"] = "large hail"
        if event_type == "wind":
            highwind = rpts.mag >= 65.
            if any(highwind):
                rpts.loc[highwind, "event_type"] = "high wind"

        # tz=3 is CST
        # tz=9 is GMT
        # To convert to UTC, add 9 and subtract tz hours.

        rpts["time"] = rpts['date_time'] + pd.to_timedelta(9 - rpts.tz, unit='h')
        # make time-zone aware datetime object
        rpts["time"] = rpts["time"].dt.tz_localize(pytz.UTC)

        if any(rpts['tz'] != 3):
            print(rpts_file, file=sys.stderr)
            unexpected_timezones = rpts[['om','date_time','tz']][rpts['tz'] != 3]
            print("get_storm_reports: found",len(unexpected_timezones),"unexpected timezones")
            if debug:
                print(unexpected_timezones, file=sys.stderr)
            print("WARNING - proceeding with program. Wrote to SPC Mar 22 2017 about fixing these lines", file=sys.stderr)

        time_window = (rpts.time >= start) & (rpts.time < end)
        rpts = rpts[time_window]
        if debug:
            print("found",len(rpts),event_type,"reports")
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
    # Use gribcode if specified. Otherwise lightning by default (187).
    # 
    # Each observation line will consist of the following 11 columns of data:
    met_columns = ["Message_Type", "Station_ID", "Valid_Time", "Lat", "Lon", "Elevation", "Grib_Code", "Level", "Height", "QC_String", "Observation_Value"]
    
    df["Message_Type"] = "ADPSFC"
    df["Station_ID"]  = df.stn
    df["Valid_Time"] = df.time.dt.strftime('%Y%m%d_%H%M%S') #df.time should be aware of UTC time zone
    df["Lat"]  = (df.slat+df.elat)/2
    df["Lon"]  = (df.slon+df.elon)/2
    df["Elevation"]  = "NA"
    df["Grib_Code"] = gribcode # 187 :lightning
    df["Level"] = "NA"
    df["Height"] = "NA"
    df["QC_String"] = "NA"
    df["Observation_Value"] = 1
    # index=False don't write index number
    return df.to_string(columns=met_columns,index=False,header=False)



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


