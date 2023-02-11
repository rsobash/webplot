import os # for NCARG_ROOT
import numpy as np
tenths = np.arange(0.1,1.1,0.1)
fifths = np.arange(0.2,1.2,0.2)

def readcm(name):
    '''Read colormap from file formatted as 0-1 RGB CSV'''
    projdir = '/glade/u/home/wrfrt/wwe/python_scripts/'
    fh = open(projdir+name, 'r')
    rgb = np.loadtxt(fh)
    fh.close()
    return rgb.tolist()

def readNCLcm(name):
    '''Read in NCL colormap for use in matplotlib
       Replaces original function in /glade/u/home/wrfrt/wwe/python_scripts/fieldinfo.py
    '''

    # comments start with ; or #
    # first real line is ncolors = 256 (or something like that)
    # The rest is bunch of rgb values, one trio per line.

    fh = open(os.getenv('NCARG_ROOT','/glade/u/apps/opt/ncl/6.5.0/intel/17.0.1')+'/lib/ncarg/colormaps/%s.rgb'%name, 'r')
    rgb = np.loadtxt(fh, comments=[';', '#', 'n']) # treat ncolors=x as a comment
    fh.close()
    if rgb.max() > 1:
        rgb = rgb/255.0
    return rgb.tolist()


fieldinfo = {
  # surface and convection-related entries
  'precip'       :{ 'levels' : [0,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.75,1,1.5,2,2.5,3.0], 'cmap': [readNCLcm('precip2_17lev')[i] for i in (0,1,2,4,5,6,7,8,10,12,13,14,15)], 'fname':['rainnc'], 'filename':'diag' },
  'precip-24hr'  :{ 'levels' : [0,0.01,0.05,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0], 'cmap': [readNCLcm('precip2_17lev')[i] for i in (0,1,2,4,5,6,7,8,10,12,13,14,15,16,17)]+['#777777', '#AAAAAA', '#CCCCCC', '#EEEEEE']+[readNCLcm('sunshine_diff_12lev')[i] for i in (4,2,1)], 'fname':['rainnc'], 'filename':'diag' },
  'precip-48hr'  :{ 'levels' : [0,0.01,0.05,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0], 'cmap': [readNCLcm('precip2_17lev')[i] for i in (0,1,2,4,5,6,7,8,10,12,13,14,15,16,17)]+['#777777', '#AAAAAA', '#CCCCCC', '#EEEEEE']+[readNCLcm('sunshine_diff_12lev')[i] for i in (4,2,1)], 'fname':['rainnc'], 'filename':'diag' },
  'precipacc'    :{ 'levels' : [0,0.01,0.05,0.1,0.25,0.5,0.75,1.0,1.25,1.5,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0], 'cmap': [readNCLcm('precip2_17lev')[i] for i in (0,1,2,4,5,6,7,8,10,12,13,14,15,16,17)]+['#777777', '#AAAAAA', '#CCCCCC', '#EEEEEE']+[readNCLcm('sunshine_diff_12lev')[i] for i in (4,2,1)], 'fname':['rainnc'], 'filename':'diag' },
  'sbcape'       :{ 'levels' : [100,250,500,750,1000,1250,1500,1750,2000,2500,3000,3500,4000,4500,5000,5500,6000],
                    'cmap'   : ['#eeeeee', '#dddddd', '#cccccc', '#aaaaaa']+readNCLcm('precip2_17lev')[3:-1], 'fname': ['sbcape'], 'filename':'diag' },
  'mlcape'       :{ 'levels' : [100,250,500,750,1000,1250,1500,1750,2000,2500,3000,3500,4000,4500,5000,5500,6000],
                    'cmap'   : ['#eeeeee', '#dddddd', '#cccccc', '#aaaaaa']+readNCLcm('precip2_17lev')[3:-1], 'fname': ['mlcape'], 'filename':'diag' },
  'mucape'       :{ 'levels' : [10,25,50,100,250,500,750,1000,1250,1500,1750,2000,2500,3000,3500,4000,4500,5000,5500,6000],
                    'cmap'   : ['#eeeeee', '#dddddd', '#cccccc', '#aaaaaa']+readNCLcm('precip2_17lev')[3:-1], 'fname': ['cape'], 'filename':'diag' },
  'sbcinh'       :{ 'levels' : [50,75,100,150,200,250,500], 'cmap': readNCLcm('topo_15lev')[1:], 'fname': ['sbcin'], 'filename':'diag' },
  'mlcinh'       :{ 'levels' : [50,75,100,150,200,250,500], 'cmap': readNCLcm('topo_15lev')[1:], 'fname': ['mlcin'], 'filename':'diag' },
  'pwat'         :{ 'levels' : [0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0],
                    'cmap'   : ['#dddddd', '#cccccc', '#e1e1d3', '#e1d5b1', '#ffffd5', '#e5ffa7', '#addd8e', '#41ab5d', '#007837', '#005529', '#0029b1'],
                    'fname'  : ['precipw'], 'filename':'diag'},
  't2'           :{ 'levels' : [-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120], 'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['t2m'], 'filename':'diag' },
  't2depart'     :{ 'levels' : [-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120], 'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['t2m'], 'filename':'diag' },
  't2-0c'         :{ 'levels' : [32], 'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['t2m'], 'filename':'diag' },
  'mslp'         :{ 'levels' : [960,964,968,972,976,980,984,988,992,996,1000,1004,1008,1012,1016,1020,1024,1028,1032,1036,1040,1044,1048,1052], 'cmap':readNCLcm('nice_gfdl')[3:193], 'fname':['mslp'], 'filename': 'diag' },
  'td2'          :{ 'levels' : [-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55,60,64,68,72,76,80,84],
                    'cmap':['#ad598a', '#c589ac','#dcb8cd','#e7cfd1','#d0a0a4','#ad5960', '#8b131d', '#8b4513','#ad7c59', '#c5a289','#dcc7b8','#eeeeee', '#dddddd', '#bbbbbb', '#e1e1d3', '#e1d5b1','#ccb77a','#ffffe5','#f7fcb9', '#addd8e', '#41ab5d', '#006837', '#004529', '#195257', '#4c787c'],
                    'fname'  : ['dewpoint_surface'], 'filename':'diag'},
  'td2depart'    :{ 'levels' : [20,25,30,35,40,45,50,55,60,64,68,72,76,80,84], 'cmap'   : ['#eeeeee', '#dddddd', '#bbbbbb', '#e1e1d3', '#e1d5b1','#ccb77a','#ffffe5','#f7fcb9', '#addd8e', '#41ab5d', '#006837', '#004529', '#195257', '#4c787c'],
                    'fname'  : ['dewpoint_surface'], 'filename':'diag'},
  'thetae'       :{  'levels' : [300,305,310,315,320,325,330,335,340,345,350,355,360], 'cmap': ['#eeeeee', '#dddddd', '#cccccc', '#aaaaaa']+readNCLcm('precip2_17lev')[3:-1], 'fname'  : ['t2m', 'q2', 'surface_pressure'], 'filename': 'diag'},
  'thetapv'      :{  'levels' : np.arange(278,386,4), 'cmap': readNCLcm('WhiteBlueGreenYellowRed'), 'fname'  : ['theta_pv'], 'filename': 'diag'},
  'rh2m'           :{  'levels' : [0,10,20,30,40,50,60,70,80,90,100,110], 'cmap': readNCLcm('precip2_17lev')[:17][::-1], 'fname': ['t2m', 'surface_pressure', 'q2'], 'filename': 'diag'},
  'pblh'         :{ 'levels' : [0,250,500,750,1000,1250,1500,1750,2000,2500,3000,3500,4000],
      'cmap': ['#eeeeee', '#dddddd', '#cccccc', '#bbbbbb', '#44aaee','#88bbff', '#aaccff', '#bbddff', '#efd6c1', '#e5c1a1', '#eebb32', '#bb9918'], 'fname': ['hpbl'], 'filename': 'diag'},
  'hmuh'         :{ 'levels' : [10,25,50,75,100,125,150,175,200,250,300,400,500], 'cmap': readNCLcm('prcp_1')[1:15], 'fname': ['updraft_helicity_max'], 'filename':'diag'},
  'hmneguh'         :{ 'levels' : [10,25,50,75,100,125,150,175,200,250,300,400,500], 'cmap': readNCLcm('prcp_1')[1:15], 'fname': ['UP_HELI_MIN'], 'filename':'diag'},
  'hmuh03'       :{ 'levels' : [10,25,50,75,100,125,150,175,200,250,300,400,500], 'cmap': readNCLcm('prcp_1')[1:15], 'fname': ['updraft_helicity_max03'], 'filename':'diag'},
  'hmuh01'       :{ 'levels' : [5,10,25,50,75,100,125,150,175,200,250,300,400,500], 'cmap': readNCLcm('prcp_1')[1:15], 'fname': ['updraft_helicity_max01'], 'filename':'diag'},
  'rvort1'       :{ 'levels' : [0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015], 'cmap': readNCLcm('prcp_1')[1:15], 'fname': ['rvort1_max'], 'filename':'diag'},
  'sspf'         :{ 'levels' : [10,25,50,75,100,125,150,175,200,250,300,400,500], 'cmap': readNCLcm('prcp_1')[:15], 'fname': ['updraft_helicity_max','WSPD10MAX','HAIL_MAXK1'], 'filename':'diag'},
  'hmup'         :{ 'levels' : [4,6,8,10,12,14,16,18,20,24,28,32,36,40,44,48], 'cmap': readNCLcm('prcp_1')[1:16], 'fname': ['w_velocity_max'], 'filename':'diag' },
  #'hmdn'         :{ 'levels' : [-19,-17,-15,-13,-11,-9,-7,-5,-3,-1,0], 'cmap': readNCLcm('prcp_1')[16:1:-1]+['#ffffff'], 'fname': ['w_velocity_min'], 'filename':'diag' },
  'hmdn'         :{ 'levels' : [2,3,4,6,8,10,12,14,16,18,20,22,24,26,28,30], 'cmap': readNCLcm('prcp_1')[1:16], 'fname': ['w_velocity_min'], 'filename':'diag' },
  #'hmwind'       :{ 'levels' : [10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42], 'cmap': readNCLcm('prcp_1')[:16], 'fname': ['WSPD10MAX'], 'filename':'diag' },
  'hmwind'       :{ 'levels' : [10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42], 'cmap': readNCLcm('prcp_1')[:16], 'fname': ['wind_speed_level1_max'], 'filename':'diag' },
  #'hmwind'       :{ 'levels' : [10,12,14,16,18,20,22,24,26,28,30,32,34], 'cmap': readNCLcm('prcp_1')[1:15], 'fname': ['WSPD10MAX'], 'filename':'diag' },
  #'hmwind'       :{ 'levels' : [20,25,30,35,40,45,50,55,60,65,70,75,80,85], 'cmap': readNCLcm('prcp_1')[1:16], 'fname': ['WSPD10MAX'], 'filename':'diag' },
  'hmgrp'        :{ 'levels' : [0.01,0.1,0.25,0.5,0.75,1.0,1.5,2.0,2.5,3.0,4.0,5.0], 'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['grpl_max'], 'filename':'diag' },
  'cref'         :{ 'levels' : [5,10,15,20,25,30,35,40,45,50,55,60,65,70], 'cmap': readcm('cmap_rad.rgb')[1:14], 'fname': ['refl10cm_max'], 'filename':'diag' },
  'ref1km'       :{ 'levels' : [5,10,15,20,25,30,35,40,45,50,55,60,65,70], 'cmap': readcm('cmap_rad.rgb')[1:14], 'fname': ['refl10cm_1km'], 'filename':'diag' },
  'srh3'         :{ 'levels' : [50,100,150,200,250,300,400,500], 'cmap': readNCLcm('perc2_9lev'), 'fname': ['srh_0_3km'], 'filename' : 'diag' },
  'srh1'         :{ 'levels' : [50,100,150,200,250,300,400,500], 'cmap': readNCLcm('perc2_9lev'), 'fname': ['srh_0_1km'], 'filename' : 'diag' },
  'shr06mag'     :{ 'levels' : [30,35,40,45,50,55,60,65,70,75,80], 'cmap': readNCLcm('perc2_9lev'), 'fname': ['uzonal_6km', 'umeridional_6km', 'uzonal_surface', 'umeridional_surface'], 'filename':'diag'},
  'shr01mag'     :{ 'levels' : [30,35,40,45,50,55,60,65,70,75,80], 'cmap': readNCLcm('perc2_9lev'), 'fname': ['uzonal_1km', 'umeridional_1km', 'uzonal_surface', 'umeridional_surface'], 'filename':'diag'},
  'zlfc'         :{ 'levels' : [0,250,500,750,1000,1250,1500,2000,2500,3000,3500,4000,5000], 'cmap': [readNCLcm('nice_gfdl')[i] for i in [3,20,37,54,72,89,106,123,141,158,175,193]], 'fname': ['lfc'], 'filename':'diag' },
  'zlcl'         :{ 'levels' : [0,250,500,750,1000,1250,1500,2000,2500,3000,3500,4000,5000], 'cmap': [readNCLcm('nice_gfdl')[i] for i in [3,20,37,54,72,89,106,123,141,158,175,193]], 'fname': ['lcl'], 'filename':'diag' },
  'ltg1'         :{ 'levels' : [0.1,0.5,1,1.5,2,2.5,3,4,5,6,7,8,10,12], 'cmap': readNCLcm('prcp_1')[:15], 'fname': ['LTG1_MAX'], 'filename':'diag' },
  'ltg2'         :{ 'levels' : [0.1,0.5,1,1.5,2,2.5,3,4,5,6,7,8,10,12], 'cmap': readNCLcm('prcp_1')[:15], 'fname': ['LTG2_MAX'], 'filename':'diag' },
  'ltg3'         :{ 'levels' : [0.1,0.5,1,1.5,2,2.5,3,4,5,6,7,8,10,12], 'cmap': readNCLcm('prcp_1')[:15], 'fname': ['LTG3_MAX'], 'filename':'diag' },
  'olrtoa'       :{ 'levels' : range(70,340,10), 'cmap': readcm('cmap_satir.rgb')[32:1:-1], 'fname': ['olrtoa'], 'filename':'diag' },

  # winter fields
  'thck1000-500' :{ 'levels' : [480,486,492,498,504,510,516,522,528,534,540,546,552,558,564,570,576,582,588,592,600], 'cmap':readNCLcm('perc2_9lev'), 'fname':['GHT_PL', 'GHT_PL'], 'arraylevel':[0,5], 'filename':'diag'}, # CSS mod
  'thck1000-850' :{ 'levels' : [82,85,88,91,94,97,100,103,106,109,112,115,118,121,124,127,130,133,136,139,142,145,148,151,154,157,160], 'cmap':readNCLcm('perc2_9lev'), 'fname':['GHT_PL', 'GHT_PL'], 'arraylevel':[0,2], 'filename':'diag'}, # CSS mod

  # pressure level entries
  'hgt200'       :{ 'levels' : list(range(10900,12500,60)),                                                                                                                    'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['height_200hPa'], 'filename':'diag'},
  'hgt250'       :{ 'levels' : [9700,9760,9820,9880,9940,10000,10060,10120,10180,10240,10300,10360,10420,10480,10540,10600,10660,10720,10780,10840,10900,10960,11020],         'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['height_250hPa'], 'filename':'diag'},
  'hgt300'       :{ 'levels' : [8400,8460,8520,8580,8640,8700,8760,8820,8880,8940,9000,9060,9120,9180,9240,9300,9360,9420,9480,9540,9600,9660,9720,9780,9840,9900,9960,10020], 'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['height_300hPa'], 'filename':'diag'},
  'hgt500'       :{ 'levels' : [4800,4860,4920,4980,5040,5100,5160,5220,5280,5340,5400,5460,5520,5580,5640,5700,5760,5820,5880,5940,6000], 'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['height_500hPa'], 'filename':'diag'},
  'hgt700'       :{ 'levels' : [2700,2730,2760,2790,2820,2850,2880,2910,2940,2970,3000,3030,3060,3090,3120,3150,3180,3210,3240,3270,3300], 'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['height_700hPa'], 'filename':'diag'},
  'hgt850'       :{ 'levels' : [1200,1230,1260,1290,1320,1350,1380,1410,1440,1470,1500,1530,1560,1590,1620,1650],                          'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['height_850hPa'], 'filename':'diag'},
  'hgt925'       :{ 'levels' : [550,580,610,640,670,700,730,760,790,820,850,880,910,940,970,1000,1030],                                    'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['height_925hPa'], 'filename':'diag'},
  # RAS: adjusted these ranges - need to capture higher wind speeds
  'speed200'     :{ 'levels' : [25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110], 'cmap': readNCLcm('wind_17lev'), 'fname': ['uzonal_200hPa','umeridional_200hPa'], 'filename':'diag'},
  'speed250'     :{ 'levels' : [25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110], 'cmap': readNCLcm('wind_17lev'), 'fname': ['uzonal_250hPa','umeridional_250hPa'], 'filename':'diag'},
  'speed300'     :{ 'levels' : [25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110], 'cmap': readNCLcm('wind_17lev'), 'fname': ['uzonal_300hPa','umeridional_300hPa'], 'filename':'diag'},
  'speed500'     :{ 'levels' : [15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100],   'cmap': readNCLcm('wind_17lev'), 'fname': ['uzonal_500hPa','umeridional_500hPa'], 'filename':'diag'},
  'speed700'     :{ 'levels' : [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85],        'cmap': readNCLcm('wind_17lev'), 'fname': ['uzonal_700hPa','umeridional_700hPa'], 'filename':'diag'},
  'speed850'     :{ 'levels' : [6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74],     'cmap': readNCLcm('wind_17lev'), 'fname': ['uzonal_850hPa','umeridional_850hPa'], 'filename':'diag'},
  'speed925'     :{ 'levels' : [6,10,14,18,22,26,30,34,38,42,46,50,54,58,62,66,70,74],     'cmap': readNCLcm('wind_17lev'), 'fname': ['uzonal_925hPa','umeridional_925hPa'], 'filename':'diag'},
  'temp200'      :{ 'levels' : [-65,-63,-61,-59,-57,-55,-53,-51,-49,-47,-45,-43,-41,-39,-37,-35,-33,-31,-29], 'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['temperature_200hPa'], 'filename':'diag'},
  'temp250'      :{ 'levels' : [-65,-63,-61,-59,-57,-55,-53,-51,-49,-47,-45,-43,-41,-39,-37,-35,-33,-31,-29], 'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['temperature_250hPa'], 'filename':'diag'},
  'temp300'      :{ 'levels' : [-65,-63,-61,-59,-57,-55,-53,-51,-49,-47,-45,-43,-41,-39,-37,-35,-33,-31,-29], 'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['temperature_300hPa'], 'filename':'diag'},
  'temp500'      :{ 'levels' : [-41,-39,-37,-35,-33,-31,-29,-26,-23,-20,-17,-14,-11,-8,-5,-2],                'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['temperature_500hPa'], 'filename':'diag'},
  'temp700'      :{ 'levels' : [-36,-33,-30,-27,-24,-21,-18,-15,-12,-9,-6,-3,0,3,6,9,12,15,18,21],            'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['temperature_700hPa'], 'filename':'diag'},
  'temp850'      :{ 'levels' : [-30,-27,-24,-21,-18,-15,-12,-9,-6,-3,0,3,6,9,12,15,18,21,24,27,30],           'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['temperature_850hPa'], 'filename':'diag'},
  'temp925'      :{ 'levels' : [-24,-21,-18,-15,-12,-9,-6,-3,0,3,6,9,12,15,18,21,24,27,30,33],                'cmap': readNCLcm('nice_gfdl')[3:193], 'fname': ['temperature_925hPa'], 'filename':'diag'},
  'td500'        :{ 'levels' : [-30,-25,-20,-15,-10,-5,0,5,10],                     'cmap' : readNCLcm('nice_gfdl')[3:193], 'fname': ['dewpoint_500hPa'], 'filename':'diag'},
  'td700'        :{ 'levels' : [-30,-25,-20,-15,-10,-5,0,5,10],                     'cmap' : readNCLcm('nice_gfdl')[3:193], 'fname': ['dewpoint_700hPa'], 'filename':'diag'},
  'td850'        :{ 'levels' : [-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30], 'cmap' : readNCLcm('nice_gfdl')[3:193], 'fname': ['dewpoint_850hPa'], 'filename':'diag'},
  'td925'        :{ 'levels' : [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30],         'cmap' : readNCLcm('nice_gfdl')[3:193], 'fname': ['dewpoint_925hPa'], 'filename':'diag'},
  'vort500'      :{ 'levels' : [0,9,12,15,18,21,24,27,30,33], 'cmap': readNCLcm('prcp_1'), 'fname': ['vorticity_500hPa'], 'filename':'diag'},
  'vort850'      :{ 'levels' : [0,9,12,15,18,21,24,27,30,33], 'cmap': readNCLcm('prcp_1'), 'fname': ['vorticity_850hPa'], 'filename':'diag'},
  'vort700'      :{ 'levels' : [0,9,12,15,18,21,24,27,30,33], 'cmap': readNCLcm('prcp_1'), 'fname': ['vorticity_700hPa'], 'filename':'diag'},
  'vortpv'       :{ 'levels' : [0,9,12,15,18,21,24,27,30,33], 'cmap': readNCLcm('prcp_1'), 'fname': ['vort_pv'], 'filename':'diag'},
  'rh300'        :{ 'levels' : [0,10,20,30,40,50,60,70,80,90,100], 'cmap' : [readNCLcm('MPL_PuOr')[i] for i in (2,18,34,50)]+[readNCLcm('MPL_Greens')[j] for j in (2,17,50,75,106,125)], 'fname': ['relhum_300hPa'], 'filename':'diag'},
  'rh500'        :{ 'levels' : [0,10,20,30,40,50,60,70,80,90,100], 'cmap' : [readNCLcm('MPL_PuOr')[i] for i in (2,18,34,50)]+[readNCLcm('MPL_Greens')[j] for j in (2,17,50,75,106,125)], 'fname': ['relhum_500hPa'], 'filename':'diag'},
  'rh700'        :{ 'levels' : [0,10,20,30,40,50,60,70,80,90,100], 'cmap' : readNCLcm('CBR_drywet'), 'fname': ['relhum_700hPa'], 'filename':'diag'},
  'rh850'        :{ 'levels' : [0,10,20,30,40,50,60,70,80,90,100], 'cmap' : readNCLcm('CBR_drywet'), 'fname': ['relhum_850hPa'], 'filename':'diag'},
  'rh925'        :{ 'levels' : [0,10,20,30,40,50,60,70,80,90,100], 'cmap' : readNCLcm('CBR_drywet'), 'fname': ['relhum_925hPa'], 'filename':'diag'},
'pvort320k'    :{ 'levels' : [0,0.1,0.2,0.3,0.4,0.5,0.75,1,1.5,2,3,4,5,7,10],
                  'cmap'   : ['#ffffff','#eeeeee','#dddddd','#cccccc','#bbbbbb','#d1c5b1','#e1d5b9','#f1ead3','#003399','#0033FF','#0099FF','#00CCFF','#8866FF','#9933FF','#660099'],
                 'fname': ['PVORT_320K'], 'filename':'upp' },
 'speed10m'     :{ 'levels' : [3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51], 'cmap': readNCLcm('wind_17lev')[1:],'fname'  : ['u10', 'v10'], 'filename':'diag'},
 'speed10m-tc'  :{ 'levels' : [6,12,18,24,30,36,42,48,54,60,66,72,78,84,90,96,102], 'cmap': readNCLcm('wind_17lev')[1:],'fname'  : ['u10', 'v10'], 'filename':'diag'},
 'stp'          :{ 'levels' : [0.5,0.75,1.0,1.5,2.0,3.0,4.0,5.0,6.0,7.0,8.0], 'cmap':readNCLcm('perc2_9lev'), 'fname':['sbcape','lcl','srh_0_1km','uzonal_6km','umeridional_6km','uzonal_surface','umeridional_surface'], 'filename':'diag'},
 'uhratio'      :{ 'levels' : [0.1,0.3,0.5,0.7,0.9,1.0,1.1,1.2,1.3,1.4,1.5], 'cmap':readNCLcm('perc2_9lev'), 'fname':['updraft_helicity_max03', 'updraft_helicity_max'], 'filename':'diag'},
 'crefuh'       :{ 'levels' : [5,10,15,20,25,30,35,40,45,50,55,60,65,70], 'cmap': readcm('cmap_rad.rgb')[0:13], 'fname': ['refl10cm_max', 'updraft_helicity_max'], 'filename':'diag' },

  # wind barb entries
  'wind10m'      :{ 'fname'  : ['u10', 'v10'], 'filename':'diag', 'skip':50 },
  'windsfc'      :{ 'fname'  : ['uzonal_surface','umeridional_surface'], 'filename': 'diag', 'skip':50 },
  'wind1km'      :{ 'fname'  : ['uzonal_1km','umeridional_1km'],         'filename': 'diag', 'skip':50 },
  'wind6km'      :{ 'fname'  : ['uzonal_6km','umeridional_6km'],         'filename': 'diag', 'skip':50 },
  'windpv'       :{ 'fname'  : ['u_pv','v_pv'],                          'filename': 'diag', 'skip':50 },
  'shr06'        :{ 'fname'  : ['uzonal_6km','umeridional_6km','uzonal_surface','umeridional_surface'], 'filename': 'diag', 'skip':50 },
  'shr01'        :{ 'fname'  : ['uzonal_1km','umeridional_1km','uzonal_surface','umeridional_surface'], 'filename': 'diag', 'skip':50 },
}
# Enter wind barb info for list of pressure levels
for plev in ['200', '250', '300', '500', '700', '850', '925']:
    fieldinfo['wind'+plev] = { 'fname' : ['uzonal_'+plev+'hPa', 'umeridional_'+plev+'hPa'], 'filename':'diag', 'skip':50}


# Combine levels from RAIN, FZRA, ICE, and SNOW for plotting 1-hr accumulated precip for each type. Ahijevych added this
#fieldinfo['ptypes']['levels'] = [fieldinfo['precip']['levels'][1:],fieldinfo['snow']['levels'],fieldinfo['ice']['levels'],fieldinfo['fzra']['levels']]

# domains = { 'domainname': { 'corners':[ll_lat,ll_lon,ur_lat,ur_lon], 'figsize':[w,h] } }
domains = { 'CONUS' :{ 'corners': [23.1593,-120.811,46.8857,-65.0212], 'fig_width': 1080 },
            'NA'  :{ 'corners': [15.00,-170.00,65.00,-50.00], 'fig_width':1080 },
            'SGP' :{ 'corners': [25.3,-107.00,36.00,-88.70], 'fig_width':1080 },
            'NGP' :{ 'corners': [40.00,-105.0,50.30,-82.00], 'fig_width':1080 },
            'CGP' :{ 'corners': [33.00,-107.50,45.00,-86.60], 'fig_width':1080 },
            'SW'  :{ 'corners': [28.00,-121.50,44.39,-102.10], 'fig_width':1080 },
            'NW'  :{ 'corners': [37.00,-124.40,51.60,-102.10], 'fig_width':1080 },
            'SE'  :{ 'corners': [26.10,-92.75,36.00,-71.00], 'fig_width':1080 },
            'NE'  :{ 'corners': [38.00,-91.00,46.80,-65.30], 'fig_width':1080 },
            'MATL':{ 'corners': [33.50,-92.25,41.50,-68.50], 'fig_width':1080 },
}
