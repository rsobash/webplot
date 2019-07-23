import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import *
from datetime import *
import pickle
import os, sys, time, argparse
import scipy.ndimage as ndimage
from scipy import interpolate
from scipy.spatial import qhull
import subprocess
import mpas, mpas_vort_cell
import re
from fieldinfo import *
# To use MFDataset, first convert netcdf4 to netcdf4-classic with nccopy
# for example, nccopy -d nc7 /glade/p/nsc/nmmm0046/schwartz/MPAS_ens_15-3km_mesh/POST/2017050100/ens_1/diag_latlon_g193.2017-05-02_00.00.00.nc /glade/scratch/ahijevyc/hwt2017/2017050100/ens_1/diag_latlon_g193.2017-05-02_00.00.00.nc
from netCDF4 import Dataset, MFDataset
# xarray can handle muliple files with NETCDF4 non-classic. netCDF4 MFDataset can't do that.
import xarray

class webPlot:
    '''A class to plot data from NCAR ensemble'''
    def __init__(self):
        self.opts = parseargs()
        self.initdate = datetime.strptime(self.opts['date'], '%Y%m%d%H')
        self.title = self.opts['title']
        self.debug = self.opts['debug']
        self.autolevels = self.opts['autolevels']
        self.domain = self.opts['domain']
        if ',' in self.opts['timerange']: self.shr, self.ehr = list(map(int, self.opts['timerange'].split(',')))
        else: self.shr, self.ehr = int(self.opts['timerange']), int(self.opts['timerange'])
        self.createFilename()
        self.ENS_SIZE = int(os.getenv('ENS_SIZE', 10))
 
    def createFilename(self):
        for f in ['fill', 'contour','barb']: # CSS added this for loop and everything in it
           if 'name' in self.opts[f]:
              if 'thresh' in self.opts[f]:
                 prefx = self.opts[f]['name']+'_'+self.opts[f]['ensprod']+'_'+str(self.opts[f]['thresh'])   # CSS
              else:
                 prefx = self.opts[f]['name']+'_'+self.opts[f]['ensprod'] # CSS
              break

        if self.shr == self.ehr:  # CSS
           self.outfile = prefx+'_f'+'%03d'%self.shr+'_'+self.domain+'.png' # 'test.png' # CSS
        else: # CSS
           #self.outfile = prefx+'_f'+'%03d'%self.shr+'-f'+'%03d'%self.ehr+'_'+self.domain+'_'+self.opts['date']+'.png' # 'test.png' # CSS
           self.outfile = prefx+'_f'+'%03d'%self.shr+'-f'+'%03d'%self.ehr+'_'+self.domain+'.png' # 'test.png' # CSS

        # create yyyymmddhh/domain/ directory if needed
        subdir_path = os.path.join(self.opts['date'], self.domain)
        if not os.path.isdir(subdir_path):
            print("webPlot.createFilename(): making new output directory "+subdir_path)
            os.mkdir(subdir_path)
        # prepend subdir_path to outfile.
        self.outfile = os.path.join(subdir_path, self.outfile)

    def toServer(self, debug=False, url="nova.mmm.ucar.edu:/web/htdocs/projects/mpas/plots/."):
        rsync_opts = '-RL'
        if debug:
            rsync_opts += 'v' # append a 'v' for verbose
        #result = subprocess.run(['rsync', rsync_opts, '--timeout=10', '--bwlimit=3', self.outfile, url], check=True, stdout=subprocess.PIPE) # tried capture_output=True but this version of subprocess doesn't recognize it.
        # send directly to koa instead of nova. Carter set up ssh keys. Not working. Get return status=14
        result = subprocess.run(['rsync', '-e',  "'ssh -vi /home/ahijevyc/.ssh/id_rsa'", '-avR', self.outfile, url], check=True, stdout=subprocess.PIPE)
        #The --contimeout option may only be used when connecting to an rsync daemon
        if result.returncode != 0:
            print(result)

    def loadMap(self, overlay=False):
        if hasattr(self, 'domain'): # once called with empty junk webplot class, just to get lats/lons/x/y/ibox attributes
            if overlay:
              PYTHON_SCRIPTS_DIR = os.getenv('PYTHON_SCRIPTS_DIR', '/glade/u/home/sobash/RT2015_gpx')   
              self.fig, self.ax, self.m = pickle.load(open('%s/overlays/rt2015_overlay_%s.pk'%(PYTHON_SCRIPTS_DIR,self.domain), 'r'))
            else:
              PYTHON_SCRIPTS_DIR = os.getenv('PYTHON_SCRIPTS_DIR', '/glade/work/ahijevyc/share/rt_ensemble/python_scripts')
              pklfile = '%s/rt2015_%s.pk'%(PYTHON_SCRIPTS_DIR,self.domain)
              print("loading "+pklfile) 
              self.fig, self.ax, self.m, self.lons,self.lats,self.min_grid_spacing_km,self.delta_deg,self.lon2d,self.lat2d,self.x2d,self.y2d,self.ibox,self.x,self.y,self.vtx,self.wts = pickle.load(open(pklfile, 'rb'))

    def readEnsemble(self):
        self.data, self.missing_members = readEnsemble(self.initdate, self.domain, timerange=[self.shr,self.ehr], fields=self.opts, debug=self.debug, ENS_SIZE=self.ENS_SIZE)

    def plotDepartures(self):
        from collections import OrderedDict
        hly_inventory, hly_forecast = OrderedDict(), OrderedDict()
        
        fieldname = self.opts['fill']['name']
        if fieldname == 't2depart': normals = open('/glade/u/home/sobash/RT2015_gpx/hly-temp-normal.txt').readlines()
        elif fieldname == 'td2depart': normals = open('/glade/u/home/sobash/RT2015_gpx/hly-dewp-normal.txt').readlines()

        # figure out local times for this UTC time because NCDC is stupid
        utcoffset = (time.mktime(time.gmtime()) - time.mktime(time.localtime()))/3600.0 #utc to local time offset for mountain time
        monthday  = (self.initdate + timedelta(hours=self.shr)).strftime(' %m %d ') #search for this string in normals file (see below)
        hr_pt        = (self.initdate + timedelta(hours=self.shr) - timedelta(hours=utcoffset+1)).hour # hour in pacific time
        hr_mt        = (self.initdate + timedelta(hours=self.shr) - timedelta(hours=utcoffset)).hour # hour in mountain time
        hr_ct        = (self.initdate + timedelta(hours=self.shr) - timedelta(hours=utcoffset-1)).hour # hour in central time
        hr_et        = (self.initdate + timedelta(hours=self.shr) - timedelta(hours=utcoffset-2)).hour # hour in eastern time 

        # get interpolated forecast at climate station locations (in F) 
        if len(self.x.shape) == 1:
            fi = interpolate.RegularGridInterpolator((self.y2d[:,0], self.x2d[0,:]), self.latlonGrid(self.data['fill'][0]), fill_value=-9999, bounds_error=False, method='linear')
        if len(self.x.shape) == 2:
            fi = interpolate.RegularGridInterpolator((self.y[:,0], self.x[0,:]), self.data['fill'][0], fill_value=-9999, bounds_error=False, method='linear')
        with open('/glade/u/home/sobash/RT2015_gpx/hly-inventory.txt') as f:
            for line in f:
                stn = line.split()
                lonlat = list(map(float, stn[1:3][::-1]))
                x_ob, y_ob = self.m(*list(zip(lonlat)))
                hly_forecast[stn[0]] = fi((y_ob,x_ob))[0]
                hly_inventory[stn[0]] = (y_ob[0], x_ob[0], lonlat[0], lonlat[1])

        # read in file of normals for this hour and plot departures
        if fieldname == 't2depart': cmap = plt.get_cmap('RdBu_r')
        elif fieldname == 'td2depart': cmap = plt.get_cmap('BrBG')
        norm = colors.BoundaryNorm([-50,-19.999,-14.999,-9.999,-4.999,0,5,10,15,20,50], cmap.N) #want to have bins so they both dont include bndrys
        self.m.set_axes_limits(ax=self.ax) # need to do this if no other basemap functions have been called
 
        for line in normals:
            if monthday in line:
                    stn = line.split() 
                    try: y_ob, x_ob, lon, lat = hly_inventory[stn[0]]
                    except: continue

                    if lon < -114: localhr = hr_pt #PT
                    elif lon > -102 and lon <= -85: localhr = hr_ct #CT
                    elif lon > -85: localhr = hr_et #ET
                    else: localhr = hr_mt #MT
                   
                    forecast  = hly_forecast[stn[0]]
                    normal    = float(stn[localhr+3][:-1])/10.0
                    departure = float(forecast - normal)
                    
                    if x_ob < self.m.xmax and x_ob > self.m.xmin and y_ob < self.m.ymax and y_ob > self.m.ymin and normal > -50:
                        if abs(departure) < 5: size = 5**2
                        elif abs(departure) < 10: size = 8**2
                        elif abs(departure) < 15: size = 11**2
                        elif abs(departure) < 20: size = 13**2
                        else: size = 15**2
                        cs = self.m.scatter(x_ob, y_ob, s=size, c=departure, cmap=cmap, norm=norm, linewidths=1.25, ax=self.ax)

        # make custom legend
        fontdict = {'family':'monospace', 'size':9 }
        x0, y0 = self.ax.transAxes.transform((0,0))
        x, y = self.fig.transFigure.inverted().transform((x0+10,y0+10))
        w, h = self.fig.get_size_inches()*self.fig.dpi
        cax = self.fig.add_axes([x,y,130/float(w),114/float(h)])
        cax.set_facecolor('#dddddd')
        cax.set_xticks([])
        cax.set_yticks([])
        for i in list(cax.spines.values()): i.set_linewidth(0.5)
 
        sizes, start_y, start_c, labels = [5,8,11,13,15], 0.84, 1, ['0-5F', '5-10F', '10-15F', '15-20F', '>= 20F']
        cax.text(0.2,0.94,'Below', va='center', ha='center', fontdict=fontdict, transform=cax.transAxes)
        cax.text(0.8,0.94,'Above', va='center', ha='center', fontdict=fontdict, transform=cax.transAxes)
        for i in range(len(sizes)):
            cax.scatter(0.2,start_y-i*0.18,s=sizes[i]**2,c=(-start_c-i*5), cmap=cmap, norm=norm, transform=cax.transAxes)
            cax.text(0.5,start_y-i*0.18,labels[i], va='center', ha='center', fontdict=fontdict, transform=cax.transAxes)
            cax.scatter(0.8,start_y-i*0.18,s=sizes[i]**2,c=(start_c+i*5), cmap=cmap, norm=norm, transform=cax.transAxes)
        
    def plotInterp(self):
        with open('stations.txt') as f: content = f.readlines()
        latlons = [ [ float(num[:-1]) for num in line.split()[-2:] ] for line in content ]
        latlons = np.array(latlons)[:,::-1]
        latlons[:,0] = -1*latlons[:,0]
        
        #latlons = zip(self.lons[10::40,10::40].flatten(), self.lats[10::40,10::40].flatten())

        f = interpolate.RegularGridInterpolator((self.y[:,0], self.x[0,:]), self.data['fill'][0], fill_value=-9999, bounds_error=False, method='linear')
        x_ob, y_ob = self.m(*list(zip(*latlons)))
        fcst_val = f((y_ob,x_ob))
 
        self.m.set_axes_limits(ax=self.ax) # need to do this if no other basemap functions have been called
        
        fontdict = {'family':'monospace', 'size':9 }
        #self.ax.cla()
        #self.ax.axis('off')
        for i in range(fcst_val.size):
            if x_ob[i] < self.m.xmax and x_ob[i] > self.m.xmin and y_ob[i] < self.m.ymax and y_ob[i] > self.m.ymin and fcst_val[i] != -9999:
                self.ax.text(x_ob[i], y_ob[i], int(round(fcst_val[i])), fontdict=fontdict, ha='center', va='center')

    def plotReports(self):
       #import csv, re, urllib2
       #url = 'http://www.spc.noaa.gov/climo/reports/%s_rpts_raw.csv'%(self.initdate.strftime('%y%m%d'))
       #print url
       #url = 'http://www.spc.noaa.gov/climo/reports/yesterday.csv'
       #url = 'http://www.spc.noaa.gov/climo/reports/today.csv'
       #response = urllib2.urlopen(url)
       #cr = csv.reader(response)
       
       import sqlite3     
       colors = ['red', 'green', 'blue']
       conn = sqlite3.connect('/glade/u/home/sobash/2013RT/REPORTS/reports_all.db')
       c = conn.cursor()      
       #for i,table in enumerate(['reports_torn', 'reports_hail', 'reports_wind']):
       for i,table in enumerate(['reports_torn', 'reports_hail', 'reports_wind']):
           if table == 'reports_hail': c.execute("SELECT slat, slon, datetime FROM %s WHERE size > 2 AND datetime > '2015-06-04 18:00' AND datetime <= '2015-06-05 00:00' ORDER BY datetime asc"%table)
           if table == 'reports_torn': c.execute("SELECT slat, slon, datetime FROM %s WHERE datetime > '2015-06-04 18:00' AND datetime <= '2015-06-05 00:00' ORDER BY datetime asc"%table)
           if table == 'reports_wind': c.execute("SELECT slat, slon, datetime FROM %s WHERE datetime > '2015-06-04 18:00' AND datetime <= '2015-06-05 00:00' ORDER BY datetime asc"%table)
           rpts = c.fetchall()
           #rpts.extend(c.fetchall())
           olats, olons, dt = zip(*rpts) 
           x_ob, y_ob = self.m(olons, olats)
           self.m.scatter(x_ob, y_ob, color=colors[i], edgecolor=None, s=25, ax=self.ax)
 
       #lats, lons, type = [], [], []
       #report_type = 'hail'
       #for row in cr:
       #    if re.search('Hail', row[0]): report_type = 'hail'
       #    if re.search('Wind', row[0]): report_type = 'wind'
       #    if re.search('Tornado', row[0]): report_type = 'torn'

       #    if re.search('^\d{4}', row[0]):
       #        lats.append(float(row[5]))
       #        lons.append(float(row[6]))
       #        type.append(report_type)
   
       #x_ob, y_ob = self.m(lons, lats)
       #self.m.scatter(x_ob, y_ob, color='k', edgecolor='k', s=3, ax=self.ax)

    def plotWarnings(self):
        try: self.m.readshapefile('/glade/u/home/sobash/SHARPpy/OBS/sbw_shp/%s'%self.initdate.strftime('%Y%m%d12'),'warnings',drawbounds=True, linewidth=1, color='black', ax=self.ax)
        except IndexError: print('IndexError - likely no warnings present')

    def plotOutlook(self):
        self.m.readshapefile('/glade/u/home/sobash/RT2015_gpx/outlook_shp/day1otlk_20150515_1200_cat','day1otlk_20150515_1200_cat',drawbounds=True, linewidth=1, color='green', ax=self.ax)

    def plotMRMS(self, thresh=40.0):
        validstr = (self.initdate+timedelta(hours=self.shr)).strftime('%Y%m%d%H')

        # READ IN MRMS GRID
        fh = Dataset('/glade/p/nmmm0001/sobash/MRMS/mrms_grid.nc', 'r')
        lats = fh.variables['lat_0'][:]
        lons = fh.variables['lon_0'][:]
        fh.close()

        # READ IN DATA
        fh = Dataset('/glade/scratch/sobash/mrms_tmp/cref_mrms_%s.nc'%validstr, 'r')
        mrms = fh.variables['CREF'][:]
        fh.close()

        lats, lons = np.meshgrid(lats, lons)
        x, y = self.m(lons, lats)

        cs1 = self.m.contourf(x, y, mrms.T, levels=[thresh,1000], colors=['k'], ax=self.ax)
        
        x0, y1 = self.ax.transAxes.transform((0,1))
        fontdict = {'family':'monospace', 'size':11}
        self.ax.text(x0+10, y1-15, 'MRMS CREF >= %d dBZ'%thresh, fontdict=fontdict, transform=None)

    def plotMRMSmax(self, field='qpe'):
        validstr = (self.initdate+timedelta(hours=self.shr-1)).strftime('%Y%m%d%H')
        #validstr = (self.initdate+timedelta(hours=self.shr-1)).strftime('%Y%m%d00')
      
        # READ IN MRMS GRID
        fh = Dataset('/glade/p/nmmm0001/sobash/MRMS/mrms_grid.nc', 'r')
        lats = fh.variables['lat_0'][:]
        lons = fh.variables['lon_0'][:]
        fh.close()

        # READ IN DATA
        if field == 'qpe':
          fname = '/glade/scratch/sobash/mrms_tmp/qpe_mrms_hrmax_%s.nc'%validstr
          field = 'QPE01'
          levels = [25.4,1000]
        if field == 'cref':
          fname = '/glade/scratch/sobash/mrms_tmp/cref_mrms_hrmax_%s.nc'%validstr 
          field = 'CREF'
          levels = [40,1000]

        fh = Dataset(fname, 'r')
        mrms = fh.variables[field][0,:]
        fh.close()

        # plot
        lats, lons = np.meshgrid(lats, lons)
        x, y = self.m(lons, lats)

        cs1 = self.m.contourf(x, y, mrms.T, levels=levels, colors=['k'], ax=self.ax)

    def plotTitleTimes(self):
        if self.opts['over']: return
        fontdict = {'family':'monospace', 'size':12, 'weight':'bold'}

        # place title and times above corners of map
        x0, y1 = self.ax.transAxes.transform((0,1))
        x0, y0 = self.ax.transAxes.transform((0,0))
        x1, y1 = self.ax.transAxes.transform((1,1))
        self.ax.text(x0, y1+10, self.title, fontdict=fontdict, transform=None)

        initstr, validstr = self.getInitValidStr() 
        self.ax.text(x1, y1+20, initstr, horizontalalignment='right', transform=None)
        self.ax.text(x1, y1+5, validstr, horizontalalignment='right', transform=None)

        # Plot missing members (use wrfout count here, if upp missing this wont show that)
        if len(self.missing_members['wrfout']) > 0:
            missing_members = sorted(set([ (x%10)+1 for x in self.missing_members['wrfout'] ])) #get member number from missing indices
            missing_members_string = ', '.join(str(x) for x in missing_members)
            self.ax.text(x1-5, y0+5, 'Missing member #s: %s'%missing_members_string, horizontalalignment='right')

    def plotFields(self):
        if 'fill' in self.data:
            if self.opts['fill']['ensprod'] == 'paintball': self.plotPaintball()
            elif self.opts['fill']['ensprod'] in ['stamp', 'maxstamp']: self.plotStamp()
            else: self.plotFill()
        
        if 'contour' in self.data:
            if self.opts['contour']['ensprod'] == 'spaghetti': self.plotSpaghetti()
            elif self.opts['contour']['ensprod'] == 'stamp': self.plotStamp()
            else: self.plotContour()
        
        if 'barb' in self.data:
            #self.plotStreamlines()
            self.plotBarbs()
       
        if self.opts['fill']['name'] in ['t2depart', 'td2depart']: self.plotDepartures()
  
    def plotFill(self, debug=False):
        if self.opts['fill']['name'] == 'ptype': self.plotFill_ptype(); return
        if self.opts['fill']['name'][0:6] == 'ptypes': self.plotFill_ptypes(); return
        elif self.opts['fill']['name'] == 'crefuh': self.plotReflectivityUH(); return

        if self.autolevels:
            min, max = self.data['fill'][0].min(), self.data['fill'][0].max()
            levels = np.linspace(min, max, num=10)
            cmap = colors.ListedColormap(self.opts['fill']['colors'])
            norm = colors.BoundaryNorm(levels, cmap.N)
            tick_labels = levels[:-1]
        else:
            levels = self.opts['fill']['levels']
            cmap = colors.ListedColormap(self.opts['fill']['colors'])
            extend, extendfrac = 'neither', 0.0
            tick_labels = levels[:-1]
            if self.opts['fill']['ensprod'] in ['prob', 'neprob', 'problt', 'probgt', 'neprobgt', 'neproblt', 'prob3d']:
                cmap = colors.ListedColormap(self.opts['fill']['colors'][:9])
                cmap.set_over(self.opts['fill']['colors'][-1])
                extend, extendfrac = 'max', 0.02
                tick_labels = levels
            norm = colors.BoundaryNorm(levels, cmap.N)

        data = self.data['fill'][0]
        # regrid 1D mesh that needs to be smoothed
        if self.opts['fill']['name'] in ['avo500', 'vort500', 'pbmin']:
            print("plotFill: regridding 1D mesh "+self.opts['fill']['name']+" to lat lon")
            data = self.latlonGrid(data)
            # smooth some of the fill fields
            if self.opts['fill']['name'] == 'avo500': data = ndimage.gaussian_filter(data, sigma=4)
            if self.opts['fill']['name'] == 'vort500': data = ndimage.gaussian_filter(data, sigma=4)
            if self.opts['fill']['name'] == 'pbmin' : data = ndimage.gaussian_filter(data, sigma=2)
            cs1 = self.m.contourf(self.x2d, self.y2d, data, levels=levels, cmap=cmap, norm=norm, extend='max', ax=self.ax)
        elif self.opts['fill']['ensprod'] in ['neprob', 'neprobgt', 'neproblt']:
            # assume the data array has been interpolated to lat/lon
            cs1 = self.m.contourf(self.x2d, self.y2d, data, levels=levels, cmap=cmap, norm=norm, extend='max', ax=self.ax)
        else:      
            # use ibox and tri
            # Sometimes you get a warning kwarg tri is ignored. 
            # Tried removing tri=True but got IndexError: too many indices for array
            #print("plotFill: starting contourf with ",self.ibox.shape," array "+self.opts['fill']['name'])
            #cs1 = self.m.contourf(self.x, self.y, self.data['fill'][0][self.ibox], tri=True, levels=levels, cmap=cmap, norm=norm, extend='max', ax=self.ax) #MPAS
            data = self.latlonGrid(self.data['fill'][0])
            if debug:
                print("plotFill: starting contourf with 2d array")
            cs1 = self.m.contourf(self.x2d, self.y2d, data, levels=levels, cmap=cmap, norm=norm, extend='max', ax=self.ax)

        self.plotColorbar(cs1, levels, tick_labels, extend, extendfrac)

    def plotFill_ptype(self):
        ml_type = np.zeros(self.data['fill'][0].shape)
        ml_type_prob = np.zeros(self.data['fill'][0].shape)
         
        for i in [1,2,3,4]:
            pts = (self.data['fill'][i-1] > ml_type_prob+0.001)
            ml_type_prob[pts] = self.data['fill'][i-1][pts]
            ml_type[pts] = i+0.001
        
        cmap = colors.ListedColormap(['#7BBF6A', 'red', 'orange', 'blue'])
        norm = colors.BoundaryNorm([1,2,3,4,5], cmap.N)
        
        x = (self.x[1:,1:] + self.x[:-1,:-1])/2.0
        y = (self.y[1:,1:] + self.y[:-1,:-1])/2.0
        cs1 = self.m.pcolormesh(x, y, np.ma.masked_equal(ml_type[1:,1:], 0), cmap=cmap, norm=norm, edgecolors='None', ax=self.ax)
        
        # make axes for colorbar, 175px to left and 30px down from bottom of map 
        x0, y0 = self.ax.transAxes.transform((0,0))
        x, y = self.fig.transFigure.inverted().transform((x0+175,y0-29.5))
        cax = self.fig.add_axes([x,y,0.985-x,y/3.0])
        cb = plt.colorbar(cs1, cax=cax, orientation='horizontal')
        cb.outline.set_linewidth(0.5)
        cb.set_ticks([0.5,1.5,2.5,3.5,4.5,5.5])
        cb.set_ticklabels(['Rain', 'Freezing Rain', 'Sleet', 'Snow'])
        cb.ax.tick_params(length=0)

    def plotFill_ptypes(self):
        types = self.data['fill']
        ntypes = len(types)

        # Plot where hourly precip of each type is greater than zero.
        # Colors match the regular ptype plot, but are translucent, so you 
        # can see them underneath each other.
        # Data arrays, colors, and levels defined in fieldinfo.py
        alpha = 0.25
        colors = self.opts['fill']['colors']
        threshes = self.opts['fill']['levels']
        type_labels= self.opts['fill']['arrayname'] # Get label directly from fieldinfo arrayname.
        type_labels = ['Rain', 'Snow', 'Sleet', 'Freezing Rain'] # Hard-coded
        if any(len(lst) != ntypes for lst in [colors, threshes, type_labels]):
            print("data, colors, threshes, and type_labels must all be same length")
            print(ntypes, colors, threshes, type_labels)
            sys.exit(1)
        # make axes for colorbar, 175px to left and 35px down from bottom of map 
        x0, y0 = self.ax.transAxes.transform((0,0))
        x, y = self.fig.transFigure.inverted().transform((x0+175,y0-35))
        # Width of space where colorbar will go. 
        cbwidth = (0.985-x)/ntypes
        for i in range(ntypes):
            levels = self.opts['fill']['levels'][i]
            # Mask pixels less than threshold 
            type = np.ma.masked_less(types[i],levels[0])
            cs = self.m.contourf(self.x, self.y, type, levels=[0,np.max(type)], colors=colors[i], ax=self.ax)
            cax = self.fig.add_axes([x+(i*cbwidth),y,0.85*cbwidth,y/3.0])
            cs.set_alpha(alpha)
            cb = plt.colorbar(cs, cax=cax, orientation='horizontal', ticks=[])
            cb.outline.set_visible(False)
            cb.ax.set_title(type_labels[i])

            data = ndimage.gaussian_filter(types[i], sigma=10)
            cs2 = self.m.contour(self.x, self.y, data, colors=colors[i], levels=levels, linewidths=1.5, ax=self.ax)
            cb = plt.colorbar(cs2, cax=cax, orientation='horizontal')
            cb.ax.tick_params(labelsize=9,length=0)
            cb.outline.set_edgecolor(colors[i])

    def plotReflectivityUH(self):
        levels = self.opts['fill']['levels']
        cmap = colors.ListedColormap(self.opts['fill']['colors'])
        norm = colors.BoundaryNorm(levels, cmap.N)
        tick_labels = levels[:-1]

        cs1 = self.m.contourf(self.x2d, self.y2d, self.latlonGrid(self.data['fill'][0]), levels=levels, cmap=cmap, norm=norm, extend='max', ax=self.ax)
        uh =  self.latlonGrid(self.data['fill'][1])
        self.m.contourf(self.x2d, self.y2d, uh, levels=[100,1000], colors='black', ax=self.ax, alpha=0.3)
        cs2 = self.m.contour( self.x2d, self.y2d, uh, levels=[100], colors='k', linewidths=0.5, ax=self.ax)
        # for some reason the zero contour is plotted if there are no other valid contours
        # are there some small negatives due to regridding? No.
        if 0.0 in cs2.levels:
            print("webplot.plotReflectivityUH has zero contour for some reason. Hide it")
            for i in cs2.collections:
                i.remove()
               
        #maxuh = self.data['fill'][1].max()
        #self.ax.text(0.03,0.03,'Domain-wide UH max %0.f'%maxuh ,ha="left",va="top",bbox=dict(boxstyle="square",lw=0.5,fc="white"), transform=self.ax.transAxes)
        
        self.plotColorbar(cs1, levels, tick_labels)

    def plotColorbar(self, cs, levels, tick_labels, extend='neither', extendfrac=0.0):
        # make axes for colorbar, 175px to left and 30px down from bottom of map 
        x0, y0 = self.ax.transAxes.transform((0,0))
        x, y = self.fig.transFigure.inverted().transform((x0+175,y0-29.5))
        cax = self.fig.add_axes([x,y,0.985-x,y/3.0])
        cb = plt.colorbar(cs, cax=cax, orientation='horizontal', extend=extend, extendfrac=extendfrac, ticks=tick_labels)
        cb.outline.set_linewidth(0.5)

    def interpolatetri(self, values, vtx, wts):
        return np.einsum('nj,nj->n', np.take(values, vtx), wts)

    def latlonGrid(self, data):
        # apply ibox to data
        data = data[self.ibox]
        if hasattr(self, "vtx") and hasattr(self, "wts"):
            data_gridded = self.interpolatetri(data, self.vtx, self.wts)
            data_gridded = np.reshape(data_gridded, self.lat2d.shape)
        else:
            print("latlonGrid: interpolating to latlon grid with griddata()")
            data_gridded = interpolate.griddata((self.lons, self.lats), data, (self.lon2d, self.lat2d), method='nearest')
        return data_gridded

    def plotContour(self):

        if self.opts['contour']['name'] in ['sbcinh','mlcinh']: linewidth, alpha = 0.5, 0.75
        else: linewidth, alpha = 1.5, 1.0
        data = self.latlonGrid(self.data['contour'][0])
        if self.opts['contour']['name'] in ['t2-0c']: data = ndimage.gaussian_filter(data, sigma=2)
        else: data = ndimage.gaussian_filter(data, sigma=25)

        cs2 = self.m.contour(self.x2d, self.y2d, data, levels=self.opts['contour']['levels'], colors='k', linewidths=linewidth, ax=self.ax, alpha=alpha)
        plt.clabel(cs2, fontsize='small', fmt='%i')

    def plotBarbs(self, debug=False):
        skip = self.opts['barb']['skip']
        if self.domain != 'CONUS': skip = int(skip*0.45)
        if self.domain == 'NA': skip = int(skip*2)

        if self.opts['fill']['name'] == 'crefuh': alpha=0.5
        else: alpha=1.0


        if debug: print("plotBarbs: starting barbs")
        # skip interval was intended for 2-D fields
        if len(self.x.shape) == 2:
            cs2 = self.m.barbs(self.x[::skip,::skip], self.y[::skip,::skip], self.data['barb'][0][::skip,::skip], self.data['barb'][1][::skip,::skip], color='black', alpha=alpha, length=5.5, linewidth=0.25, sizes={'emptybarb':0.05}, ax=self.ax)
        if len(self.x.shape) == 1:
            u2d = self.latlonGrid(self.data['barb'][0])
            v2d = self.latlonGrid(self.data['barb'][1])
            # rotate vectors so they represent the direction properly on the map projection
            u10_rot, v10_rot, x, y = self.m.rotate_vector(u2d, v2d, self.lon2d, self.lat2d, returnxy=True)
            #cs2 = self.m.barbs(self.x2d[::skip,::skip], self.y2d[::skip,::skip], u2d[::skip,::skip], v2d[::skip,::skip], color='black', alpha=alpha, length=5.5, linewidth=0.25, sizes={'emptybarb':0.05}, ax=self.ax)
            cs2 = self.m.barbs(x[::skip,::skip], y[::skip,::skip], u10_rot[::skip,::skip], v10_rot[::skip,::skip], color='black', alpha=alpha, length=5.5, linewidth=0.25, sizes={'emptybarb':0.05}, ax=self.ax)
    
    def plotStreamlines(self):
        speed = np.sqrt(self.data['barb'][0]**2 + self.data['barb'][1]**2)
        lw = 5*speed/speed.max()
        cs2 = self.m.streamplot(self.x[0,:], self.y[:,0], self.data['barb'][0], self.data['barb'][1], color='k', density=3, linewidth=lw, ax=self.ax)
        cs2.lines.set_alpha(0.5)
        cs2.arrows.set_alpha(0.5) #apparently this doesn't work?

    def plotPaintball(self):
       rects, labels = [], []
       colorlist = self.opts['fill']['colors']
       levels = self.opts['fill']['levels']
       for i in range(self.data['fill'][0].shape[0]):
           cs = self.m.contourf(self.x, self.y, self.data['fill'][0][i,self.ibox], tri=True, levels=levels, colors=[colorlist[i%len(colorlist)]], ax=self.ax, alpha=0.5)
           rects.append(plt.Rectangle((0,0),1,1,fc=colorlist[i%len(colorlist)]))
           labels.append("member %d"%(i+1))

       plt.legend(rects, labels, ncol=5, loc='right', bbox_to_anchor=(1.0,-0.05), fontsize=11, \
                  frameon=False, borderpad=0.25, borderaxespad=0.25, handletextpad=0.2)

    def plotSpaghetti(self):
       proxy = []
       colorlist = self.opts['contour']['colors']
       levels = self.opts['contour']['levels']
       data = ndimage.gaussian_filter(self.data['contour'][0], sigma=[0,4,4])
       for i in range(self.data['contour'][0].shape[0]):
           #cs = self.m.contour(self.x, self.y, data[i,:], levels=levels, colors=[colorlist[i]], linewidths=2, linestyles='solid', ax=self.ax)
           cs = self.m.contour(self.x, self.y, data[i,:], levels=levels, colors='k', alpha=0.6, linewidths=1, linestyles='solid', ax=self.ax)
           #proxy.append(plt.Rectangle((0,0),1,1,fc=colorlist[i]))
       #plt.legend(proxy, ["member %d"%i for i in range(1,11)], ncol=5, loc='right', bbox_to_anchor=(1.0,-0.05), fontsize=11, \
       #           frameon=False, borderpad=0.25, borderaxespad=0.25, handletextpad=0.2)

    def plotStamp(self, debug=False):
       fig_width_px, dpi = 1280, 90
       fig = plt.figure(dpi=dpi)

       num_rows, num_columns = 3, 4
       fig_width = fig_width_px/dpi
       width_per_panel = fig_width/float(num_columns)
       height_per_panel = width_per_panel*self.m.aspect
       fig_height = height_per_panel*num_rows
       fig_height_px = fig_height*dpi
       fig.set_size_inches((fig_width, fig_height))

       levels = self.opts['fill']['levels']
       cmap = colors.ListedColormap(self.opts['fill']['colors'])
       norm = colors.BoundaryNorm(levels, cmap.N)
       filename = self.opts['fill']['filename']

       memberidx = 0 
       for j in range(0,num_rows):
           for i in range(0,num_columns):
               member = num_columns*j+i
               if member > 9: break
               spacing_w, spacing_h = 5/float(fig_width_px), 5/float(fig_height_px)
               spacing_w = 10/float(fig_width_px)
               x, y = i*width_per_panel/float(fig_width), 1.0 - (j+1)*height_per_panel/float(fig_height)
               w, h = (width_per_panel/float(fig_width))-spacing_w, (height_per_panel/float(fig_height))-spacing_h
               if member == 9: y = 0

               if debug:
                   print('member', member, 'creating axes at', x, y)
               thisax = fig.add_axes([x,y,w,h])

               thisax.axis('on')
               for axis in ['top','bottom','left','right']: thisax.spines[axis].set_linewidth(0.5)
               self.m.drawcoastlines(ax=thisax, linewidth=0.3)
               self.m.drawstates(linewidth=0.15, ax=thisax)
               self.m.drawcountries(ax=thisax, linewidth=0.3)
               thisax.text(0.03,0.97,member+1,ha="left",va="top",bbox=dict(boxstyle="square",lw=0.5,fc="white"), transform=thisax.transAxes)
               
               # plot, unless file that has fill field is missing, then skip
               if member not in self.missing_members[filename] and member < self.ENS_SIZE:
                   data = self.latlonGrid(self.data['fill'][0][memberidx,:])
                   if debug:
                       print("plotStamp: starting contourf with regridded array", memberidx)
                   cs1 = self.m.contourf(self.x2d, self.y2d, data, levels=levels, cmap=cmap, norm=norm, extend='max', ax=thisax)
                   memberidx += 1

       # use every other tick for large colortables, remove last tick label for both
       if self.opts['fill']['name'] in ['goesch3', 'goesch4', 't2', 'precipacc' ]: ticks = levels[:-1][::2] # CSS added precipacc       
       else: ticks = levels[:-1] 

       # add colorbar to figure
       cax = fig.add_axes([0.51,0.3,0.48,0.02])
       cb = plt.colorbar(cs1, cax=cax, orientation='horizontal', ticks=ticks, extendfrac=0.0)
       cb.outline.set_linewidth(0.5)
       cb.ax.tick_params(labelsize=9)

       # add init/valid text
       fontdict = {'family':'monospace', 'size':13, 'weight':'bold'}
       initstr, validstr = self.getInitValidStr() 

       fig.text(0.51, 0.22, self.title, fontdict=fontdict, transform=fig.transFigure)
       fig.text(0.51, 0.22 - 25/float(fig_height_px), initstr, transform=fig.transFigure)
       fig.text(0.51, 0.22 - 40/float(fig_height_px), validstr, transform=fig.transFigure)

       # add NCAR logo and text below logo
       x, y = fig.transFigure.transform((0.51,0))
       fig.figimage(plt.imread('ncar.png'), xo=x, yo=y+15, zorder=1000)
       plt.text(x+10, y+5, 'ensemble.ucar.edu', fontdict={'size':9, 'color':'#505050'}, transform=None)


    def getInitValidStr(self):
       initstr  = self.initdate.strftime(' Init: %a %Y-%m-%d %H UTC')
       if ((self.ehr - self.shr) == 0):
            validstr = (self.initdate+timedelta(hours=self.shr)).strftime('Valid: %a %Y-%m-%d %H UTC')
       else:
            # match precip or precip-24hr, but not precipacc
            # accept precip-24hr, precip-48hr, precip-120hr, etc.
            if self.opts['fill']['name'] == 'precip' or is_precip_diff(self.opts['fill']['name']):
                # do not subtract 1 from start hour if array is difference of accumulated precipitation
                validstr1 = (self.initdate+timedelta(hours=(self.shr))).strftime('%a %Y-%m-%d %H UTC')
            else:
                validstr1 = (self.initdate+timedelta(hours=(self.shr-1))).strftime('%a %Y-%m-%d %H UTC')
            validstr2 = (self.initdate+timedelta(hours=self.ehr)).strftime('%a %Y-%m-%d %H UTC')
            validstr = "Valid: %s - %s"%(validstr1, validstr2)
       return initstr, validstr

    def saveFigure(self, trans=False):
        # place NCAR logo 57 pixels below bottom of map, then save image 
        if 'ensprod' in self.opts['fill']:  # CSS needed incase not a fill plot
           if not trans and self.opts['fill']['ensprod'] not in ['stamp', 'maxstamp']:
             x, y = self.ax.transAxes.transform((0,0))
             self.fig.figimage(plt.imread('ncar.png'), xo=x, yo=(y-44), zorder=1000)
             plt.text(x+10, y-54, 'ensemble.ucar.edu', fontdict={'size':9, 'color':'#505050'}, transform=None)

        plt.savefig(self.outfile, dpi=90, transparent=trans)
        
        if self.opts['convert']:
            #command = 'convert -colors 255 %s %s'%(self.outfile, self.outfile)
            if not self.opts['fill']: ncolors = 48 #if no fill field exists
            elif self.opts['fill']['ensprod'] in ['prob', 'neprob', 'probgt', 'problt', 'neprobgt', 'neproblt']: ncolors = 48
            elif self.opts['fill']['name'] in ['crefuh']: ncolors = 48
            else: ncolors = 255
            #command = '/glade/u/home/sobash/pngquant/pngquant %d %s --ext=.png --force'%(ncolors,self.outfile)
            if os.environ['NCAR_HOST'] == "cheyenne":
                bindir= '/glade/u/home/ahijevyc/bin_cheyenne/'
            else:
                bindir= '/glade/u/home/ahijevyc/bin/'
            command = bindir + 'pngquant %d %s --ext=.png --force'%(ncolors,self.outfile)
            ret = subprocess.check_call(command.split())
        plt.clf()

def parseargs():
    '''Parse arguments and return dictionary of fill, contour and barb field parameters'''

    parser = argparse.ArgumentParser(description='Web plotting script for NCAR ensemble')
    parser.add_argument('-d', '--date', default=datetime.utcnow().strftime('%Y%m%d00'), help='initialization datetime (YYYYMMDDHH)')
    parser.add_argument('-tr', '--timerange', required=True, help='time range of forecasts (START,END)')
    parser.add_argument('-f', '--fill', help='fill field (FIELD_PRODUCT_THRESH), field keys:'+','.join(list(fieldinfo.keys())))
    parser.add_argument('-c', '--contour', help='contour field (FIELD_PRODUCT_THRESH)')
    parser.add_argument('-b', '--barb', help='barb field (FIELD_PRODUCT_THRESH)')
    parser.add_argument('-bs', '--barbskip', help='barb skip interval')
    parser.add_argument('-t', '--title', help='title for plot')
    parser.add_argument('-dom', '--domain', default='CONUS', help='domain to plot')
    parser.add_argument('-al', '--autolevels', action='store_true', help='use min/max to determine levels for plot')
    parser.add_argument('-con', '--convert', default=True, action='store_false', help='run final image through imagemagick')
    parser.add_argument('-i', '--interp', default=False, action='store_true', help='plot interpolated station values')
    parser.add_argument('-sig', '--sigma', default=2, help='smooth probabilities using gaussian smoother')
    parser.add_argument('--debug', action='store_true', help='turn on debugging')
    parser.add_argument('-v', '--verif', default=None, help='plot verification data')
    parser.add_argument('--over', default=False, action='store_true', help='plot as overlay (no lines, transparent, no convert)')

    opts = vars(parser.parse_args())
    if opts['interp']: opts['over'] = True
    
    # opts = { 'date':date, 'timerange':timerange, 'fill':'sbcape_prob_25', 'ensprod':'mean' ... }
    # now, convert underscore delimited fill, contour, and barb args into dicts
    for f in ['contour','barb','fill']:
        thisdict = {}
        if opts[f] is not None:
            input = opts[f].lower().split('_')

            thisdict['name']      = input[0]
            thisdict['ensprod']   = input[1]
            thisdict['arrayname'] = fieldinfo[input[0]]['fname']
            
            # assign contour levels and colors
            if (input[1] in ['prob', 'neprob', 'probgt', 'problt', 'neprobgt', 'neproblt', 'prob3d']):
                if len(input)<3:
                    print("your -f option has less than 3 components. It needs name, ensprod, and thresh.")
                    sys.exit(1)
                thisdict['thresh']  = float(input[2])
                if int(opts['sigma']) != 40: thisdict['levels']  = np.arange(0.1,1.1,0.1)
                else: thisdict['levels']  = [0.02,0.05,0.1,0.15,0.2,0.25,0.35,0.45,0.6]
                #thisdict['levels'] = np.arange(0.1,1.1,0.2)
                thisdict['colors']  = readNCLcm('perc2_9lev')
                #thisdict['colors']  = ['#d9d9d9', '#bdbdbd', '#969696', '#636363', '#252525'] # greyscale
            elif (input[1] in ['paintball', 'spaghetti']):
                thisdict['thresh']  = float(input[2])
                thisdict['levels']  = [float(input[2]), 1000]
                thisdict['colors']  = readNCLcm('GMT_paired') 
            elif (input[1] == 'var'):
                if (input[0][0:3] == 'hgt'):
                    thisdict['levels']  = [2,4,6,8,10,15,20,25,30,35,40,45,50,55,60,65,70,75] #hgt 
                    thisdict['colors']  = readNCLcm('wind_17lev')
                elif (input[0][0:3] == 'spe'):
                    thisdict['levels']  = [1,2,3,4,5,6,7,8,9,10,12.5,15,20,25,30,35,40,45] #iso
                    thisdict['colors']  = readNCLcm('wind_17lev')
                else:
                    thisdict['levels']  = [0.5,1,1.5,2,3,4,5,6,7,8,10] #tmp/td 
                    thisdict['colors']  = readNCLcm('perc2_9lev')
            elif 'levels' in fieldinfo[input[0]]:
                thisdict['levels']  = fieldinfo[input[0]]['levels']
                thisdict['colors']  = fieldinfo[input[0]]['cmap']
          
            # get vertical array index for 3D array fields
            if 'arraylevel' in fieldinfo[input[0]]:
                thisdict['arraylevel'] = fieldinfo[input[0]]['arraylevel']
            
            # get barb-skip for barb fields
            if opts['barbskip'] is not None:    thisdict['skip'] = int(opts['barbskip'])
            elif 'skip' in fieldinfo[input[0]]: thisdict['skip'] = fieldinfo[input[0]]['skip']
            
            # get filename
            if 'filename' in fieldinfo[input[0]]: thisdict['filename'] = fieldinfo[input[0]]['filename']
            else:                                 thisdict['filename'] = 'wrfout'

        opts[f] = thisdict
    return opts

def makeEnsembleList(wrfinit, timerange, ENS_SIZE):
    # create lists of files (and missing file indices) for various file types
    shr, ehr = timerange
    file_list    = { 'wrfout':[], 'upp': [], 'diag':[] }
    missing_list = { 'wrfout':[], 'upp': [], 'diag':[] }

    EXP_DIR = os.getenv('EXP_DIR', '/glade/scratch/wrfrt/realtime_ensemble/ensf')

    missing_index = 0
    for hr in range(shr,ehr+1):
            wrfvalidstr = (wrfinit + timedelta(hours=hr)).strftime('%Y-%m-%d_%H:%M:%S')
            yyyymmddhh = wrfinit.strftime('%Y%m%d%H')
            for mem in range(1,ENS_SIZE+1):
                wrfout = '%s/%s/wrf_rundir/ens_%d/wrfout_d02_%s'%(EXP_DIR,yyyymmddhh,mem,wrfvalidstr)
                diag   = '%s/%s/wrf_rundir/ens_%d/diags_d02.%s.nc'%(EXP_DIR,yyyymmddhh,mem,wrfvalidstr)
                #diag   = '/glade/scratch/sobash/FOR_MORRIS/%s/mem%d/diags_d02_f%03d.nc'%(yyyymmddhh,mem,hr)
                diag   = '/glade/scratch/sobash/FOR_MORRIS/%s/mem%d/diags_d02_f%03d.nc'%(yyyymmddhh,mem,hr)
                upp    = '%s/%s/post_rundir/mem_%d/fhr_%d/WRFTWO%02d.nc'%(EXP_DIR,yyyymmddhh,mem,hr,hr)
                #ens1   = '/glade/p/nmm0001/romine/rt2015/ens_1km/%s/mem%02d_%s.nc'%(yyyymmddhh,mem,yyyymmddhh)
                if os.path.exists(wrfout): file_list['wrfout'].append(wrfout)
                else: missing_list['wrfout'].append(missing_index)
                if os.path.exists(diag): file_list['diag'].append(diag)
                else: missing_list['diag'].append(missing_index)
                if os.path.exists(upp): file_list['upp'].append(upp)
                else: missing_list['upp'].append(missing_index)
                missing_index += 1
    return (file_list, missing_list)

def makeEnsembleListStan(wrfinit, timerange, ENS_SIZE):
    # create lists of files (and missing file indices) for various file types
    shr, ehr = timerange
    file_list    = { 'wrfout':[], 'upp': [], 'diag':[] }
    missing_list = { 'wrfout':[], 'upp': [], 'diag':[] }

    EXP_DIR = os.getenv('EXP_DIR', '/glade/scratch/trier/jun4-5/')

    missing_index = 0
    #for hr in range(shr,ehr+1):
    print(shr, ehr)
    for m in range(shr*60,(ehr*60)+1,15):
            #wrfvalidstr = (wrfinit + timedelta(hours=m)).strftime('%Y-%m-%d_%H:%M:%S')
            wrfvalidstr = (wrfinit + timedelta(minutes=m)).strftime('%Y-%m-%d_%H:%M:%S')
            yyyymmddhh = wrfinit.strftime('%Y%m%d%H')
            for mem in range(1,ENS_SIZE+1):
                wrfout = '%s/ens_%d/wrfout_d02_%s'%(EXP_DIR,mem,wrfvalidstr)
                diag   = '%s/ens_%d/diags_d02.%s.nc'%(EXP_DIR,mem,wrfvalidstr)
                print(diag)
                if os.path.exists(wrfout): file_list['wrfout'].append(wrfout)
                else: missing_list['wrfout'].append(missing_index)
                if os.path.exists(diag): file_list['diag'].append(diag)
                else: missing_list['diag'].append(missing_index)
                missing_index += 1
    return (file_list, missing_list)

def makeEnsembleListHREF(wrfinit, timerange, ENS_SIZE):
    # create lists of files (and missing file indices) for various file types
    shr, ehr = timerange
    file_list    = { 'wrfout':[], 'diag':[] }
    missing_list = { 'wrfout':[], 'diag':[] }
    missing_index = 0
    wrfinit_prev = (wrfinit - timedelta(hours=12))

    for hr in range(shr,ehr+1):
        yyyymmddhh = wrfinit.strftime('%Y%m%d%H')
        yyyymmddhh_p = wrfinit_prev.strftime('%Y%m%d%H')
        hr_p = hr - 12
        init = wrfinit.strftime('%H')
        init_p = wrfinit_prev.strftime('%H')

        for mem in range(1,ENS_SIZE+1):
            if mem == 1: diag   = '/glade/scratch/sobash/HREF/%s/hiresw.t%sz.arw_5km.f%02d.conus.grib2'%(yyyymmddhh,init,hr)
            if mem == 2: diag   = '/glade/scratch/sobash/HREF/%s/hiresw.t%sz.arw_5km.f%02d.conusmem2.grib2'%(yyyymmddhh,init,hr)
            if mem == 3: diag   = '/glade/scratch/sobash/HREF/%s/hiresw.t%sz.nmmb_5km.f%02d.conus.grib2'%(yyyymmddhh,init,hr)
            if mem == 4: diag   = '/glade/scratch/sobash/HREF/%s/nam.t%sz.conusnest.hiresf%02d.tm00.grib2'%(yyyymmddhh,init,hr)
            if mem == 5: diag   = '/glade/scratch/sobash/HREF/%s/hiresw.t%sz.arw_5km.f%02d.conus.grib2'%(yyyymmddhh_p,init_p,hr_p)
            if mem == 6: diag   = '/glade/scratch/sobash/HREF/%s/hiresw.t%sz.arw_5km.f%02d.conusmem2.grib2'%(yyyymmddhh_p,init_p,hr_p)
            if mem == 7: diag   = '/glade/scratch/sobash/HREF/%s/hiresw.t%sz.nmmb_5km.f%02d.conus.grib2'%(yyyymmddhh_p,init_p,hr_p)
            if mem == 8: diag   = '/glade/scratch/sobash/HREF/%s/nam.t%sz.conusnest.hiresf%02d.tm00.grib2'%(yyyymmddhh_p,init_p,hr_p)

            print(diag)

            if os.path.exists(diag): file_list['diag'].append(diag)
            else: missing_list['diag'].append(missing_index)
            missing_index += 1


    return (file_list, missing_list)

def makeEnsembleListMPAS(wrfinit, timerange, ENS_SIZE, g193=False, debug=False):
    # create lists of files (and missing file indices) for various file types
    shr, ehr = timerange
    file_list    = { 'wrfout':[], 'diag':[] }
    missing_list = { 'wrfout':[], 'diag':[] }
    missing_index = 0
    for hr in range(shr,ehr+1):
        wrfvalidstr = (wrfinit + timedelta(hours=hr)).strftime('%Y-%m-%d_%H.%M.%S')
        yyyymmddhh = wrfinit.strftime('%Y%m%d%H')
        for mem in range(1,ENS_SIZE+1):
            diag   = '/glade/p/nsc/nmmm0046/schwartz/MPAS_ens_15-3km_mesh/POST/%s/ens_%d/diag.%s.nc'%(yyyymmddhh,mem,wrfvalidstr)
            if g193:
                diag   = '/glade/p/nsc/nmmm0046/schwartz/MPAS_ens_15-3km_mesh/POST/%s/ens_%d/diag_latlon_g193.%s.nc'%(yyyymmddhh,mem,wrfvalidstr)
            if debug: print(diag)
            if os.path.exists(diag): file_list['diag'].append(diag)
            else: missing_list['diag'].append(missing_index)
            missing_index += 1
    if not file_list['diag']:
        print('Empty file_list')
    return (file_list, missing_list)

def makeEnsembleListArchive(wrfinit, timerange):
    # create lists of files (and missing file indices) for various file types
    shr, ehr = timerange
    file_list    = { 'wrfout':[], 'upp': [], 'diag':[] }
    missing_list = { 'wrfout':[], 'upp': [], 'diag':[] }

    missing_index = 0
    yyyymmddhh = wrfinit.strftime('%Y%m%d%H')
    for mem in range(1,11):
        ens1   = '/glade/scratch/sobash/RT2015/%s/mem%d_surrogate_%s.nc'%(yyyymmddhh,mem,yyyymmddhh)
        #ens1   = '/glade/scratch/sobash/RT2013_1KMENS/%s/mem%02d_%s.nc'%(yyyymmddhh,mem,yyyymmddhh)
        print(ens1)
        if os.path.exists(ens1): file_list['wrfout'].append(ens1)
        else: missing_list['wrfout'].append(missing_index)
        if os.path.exists(ens1): file_list['diag'].append(ens1)
        else: missing_list['diag'].append(missing_index)
        if os.path.exists(ens1): file_list['upp'].append(ens1)
        else: missing_list['upp'].append(missing_index)
        missing_index += 1
    return (file_list, missing_list)

def makeEnsembleListNSC(wrfinit, timerange):
    # create lists of files (and missing file indices) for various file types
    shr, ehr = timerange
    file_list    = { 'wrfout':[], 'diag':[] }
    missing_list = { 'wrfout':[], 'diag':[] }
    
    EXP_DIR = os.getenv('EXP_DIR', '/glade/scratch/wrfrt/realtime_ensemble/ensf')

    missing_index = 0
    for hr in range(shr,ehr+1):
            wrfvalidstr = (wrfinit + timedelta(hours=hr)).strftime('%Y-%m-%d_%H_%M_%S')
            yyyymmddhh = wrfinit.strftime('%Y%m%d%H')
            for mem in range(1,2):
                diag   = '%s/%s/diags_d01_%s.nc'%(EXP_DIR,yyyymmddhh,wrfvalidstr)
                print(diag)
                if os.path.exists(diag): file_list['diag'].append(diag)
                else: missing_list['diag'].append(missing_index)
                missing_index += 1
    return (file_list, missing_list)

def makeEnsembleListDA(wrfinit, timerange):
    # create lists of files (and missing file indices) for various file types
    shr, ehr = timerange
    file_list    = { 'wrfout':[], 'diag':[] }
    missing_list = { 'wrfout':[], 'diag':[] }

    #EXP_DIR = os.getenv('EXP_DIR', '/glade/scratch/hclin/CONUS/wrfda/expdir/rt/fcst_15km')
    #EXP_DIR = os.getenv('EXP_DIR', '/glade/scratch/sobash/VSE/1km_pbl7')
    EXP_DIR = os.getenv('EXP_DIR', '/glade/scratch/schwartz/VSE/3km_pbl7')
    missing_index = 0
    for hr in range(shr,ehr+1):
            wrfvalidstr = (wrfinit + timedelta(hours=hr)).strftime('%Y-%m-%d_%H:%M:%S')
            yyyymmddhh = wrfinit.strftime('%Y%m%d%H')
            for mem in range(1,2):
                wrfout = '%s/%s/wrfout_d01_%s'%(EXP_DIR,yyyymmddhh,wrfvalidstr)
                #diag   = '%s/%s/wrf/join/vse_d01.%s.nc'%(EXP_DIR,yyyymmddhh,wrfvalidstr)
                diag   = '%s/%s/wrf/diags_d01.%s.nc'%(EXP_DIR,yyyymmddhh,wrfvalidstr)
                print(diag)
                if os.path.exists(wrfout): file_list['wrfout'].append(wrfout)
                else: missing_list['wrfout'].append(missing_index)
                if os.path.exists(diag): file_list['diag'].append(diag)
                else: missing_list['diag'].append(missing_index)
                missing_index += 1
    return (file_list, missing_list)

def readEnsemble(wrfinit, domain, timerange=None, fields=None, debug=False, ENS_SIZE=10):
    ''' Reads in desired fields and returns 2-D arrays of data for each field (barb/contour/field) '''
    if debug: 
        print(fields)
    
    datadict = {}
    #file_list, missing_list = makeEnsembleList(wrfinit, timerange, ENS_SIZE) #construct list of files
    #file_list, missing_list = makeEnsembleListNSC(wrfinit, timerange) #construct list of files
    #file_list, missing_list = makeEnsembleListStan(wrfinit, timerange, ENS_SIZE) #construct list of files
    file_list, missing_list = makeEnsembleListMPAS(wrfinit, timerange, ENS_SIZE, g193=False, debug=debug) #construct list of files
    #file_list, missing_list = makeEnsembleListArchive(wrfinit, timerange) #construct list of files
    #file_list, missing_list = makeEnsembleListHybrid(wrfinit, timerange) #construct list of files
    #file_list, missing_list = makeEnsembleListHREF(wrfinit, timerange, ENS_SIZE) #construct list of files
 
    # loop through fill field, contour field, barb field and retrieve required data
    for f in ['fill', 'contour', 'barb']:
        if not list(fields[f].keys()): continue
        if debug: print('Reading field:', fields[f]['name'], 'from', fields[f]['filename'])
        
        # save some variables for use in this function
        filename = fields[f]['filename']
        arrays = fields[f]['arrayname']
        fieldtype = fields[f]['ensprod']
        fieldname = fields[f]['name']
        if fieldtype in ['prob', 'neprob', 'problt', 'probgt', 'neprobgt', 'neproblt', 'prob3d']: thresh = fields[f]['thresh']
        if fieldtype[0:3]=='mem': member = int(fieldtype[3:])
        
        # open Multi-file netcdf dataset
        if debug:
            print("opening xarray mfdataset ", file_list[filename])

        fh = xarray.open_mfdataset(file_list[filename],concat_dim='Time')
        # This concatenation dimension includes different times AND members.
        fh = fh.rename({'Time':'TimeMember'})
       
        # loop through each field, wind fields will have two fields that need to be read
        datalist = []
        for n,array in enumerate(arrays):
            if debug: print('Reading', array)

            #read in 3D array (times*members,ny,nx) from file object
            if 'arraylevel' in fields[f]:
                if isinstance(fields[f]['arraylevel'], list): level = fields[f]['arraylevel'][n]
                else: level = fields[f]['arraylevel']
            else: level = None
            
            #if level == 'max':  data = np.amax(fh.variables[array][:,:,:,:], axis=1)
            #elif level is None: data = fh.variables[array][:,:,:]
            #else:               data = fh.variables[array][:,level,:,:]

            data = fh.variables[array][:,:]
            data = data.values # use the numpy array, not the full xarray object.
            # Many things that come afterward assume a numpy array, like flatten method. 
            if fh.variables[array].dims[1] == 'nVertices':
                if debug:
                    print("field on vertices, like vorticity_500hPa, put on cells")
                fieldv = data
                nEdgesOnCell, verticesOnCell = readMPASVertices()
                # verticesOnCell is the transpose of what mpas_vort_cell1 expects
                fieldc = mpas_vort_cell.mpas_vort_cell1(nEdgesOnCell, verticesOnCell.T, fieldv)
                data = fieldc 


            #data = data.reshape((10,data.shape[1]))
            #data = np.swapaxes(data,0,1) # flip first two axes so time is first
            
            # if all times are in one file, then need to reshape and extract desired times
            #data = data.reshape((10,49,*spatial_dimensions)) # reshape
            #data = data[:,timerange[0]:timerange[1]+1,:,:] # extract desired times
            #data = np.swapaxes(data,0,1) # flip first two axes so time is first
            #data = data.reshape((10*((timerange[1]+1)-timerange[0])),data.shape[2],data.shape[3]) #reshape again

            

            # change units for certain fields
            if array in ['u10','v10']:  data = data*1.93 # m/s > kt
            elif re.compile("^uzonal_\d\d+[A-Za-z]+").match(array):      data = data*1.93 # m/s > kt
            elif re.compile("^umeridional_\d\d+[A-Za-z]+").match(array): data = data*1.93 # m/s > kt
            elif re.compile("^[uv]_pv$").match(array):                   data = data*1.93 # m/s > kt
            elif array in ['dewpoint_surface', 't2m']:             data = (data - 273.15)*1.8 + 32.0 # K > F 
            elif array in ['precipw']:                             data = data*0.0393701 # mm > in             
            elif array in ['rainnc', 'grpl_max', 'SNOW_ACC_NC']:   data = data*0.0393701 # mm > in 
            elif array in ['T_PL', 'TD_PL', 'SFC_LI']:             data = data - 273.15 # K > C
            elif re.compile("^temperature_\d\d+hPa$").match(array):data = data - 273.15 # K > C
            elif re.compile("^dewpoint_\d\d+hPa$").match(array):   data = data - 273.15 # K > C
            elif array in ['mslp']:                                data = data*0.01 # Pa > hPa
            elif array in ['UP_HELI_MIN']:                         data = np.abs(data)
            elif array in ['w_velocity_min']:                      data = data*-1.0 
            elif array in ['PVORT_320K']:                          data = data*1000000 # multiply by 1e6
            elif re.compile("^vorticity_\d\d+hPa").match(array):   data = data * 1e5
            elif array == "vort_pv":                               data = data * 1e5
            elif array in ['PBMIN', 'PBMIN_SFC', 'BESTPBMIN', 'MLPBMIN', 'MUPBMIN']:  data = data*0.01 #  Pa -> hPa
#            elif array in ['LTG1_MAX1', 'LTG2_MAX', 'LTG3_MAX']:   data = data*0.20 #  scale down excess values
            
            datalist.append(data)

        # these are derived fields, we don't have in any of the input files but we can compute
        print(datalist[0].shape)
        if 'name' in fields[f]:
            if fieldname in ['shr06mag', 'shr01mag']:
                # derive wind shear from top and bottom level
                # Assume datalist is 4-element list: [bottom_u, bottom_v, top_u, top_v]
                datalist = [np.sqrt((datalist[0]-datalist[2])**2 + (datalist[1]-datalist[3])**2)]
            elif fieldname[0:5] == 'speed' or fieldname in ['bunkmag']:
                # derive speed from u and v components. different than wrf. It already has speed in S_PL array.
                datalist = [np.sqrt(datalist[0]**2 + datalist[1]**2)]
            elif fieldname in ['shr06','shr01']: datalist = [datalist[0]-datalist[2], datalist[1]-datalist[3]]
            elif fieldname == 'uhratio': datalist = [compute_uhratio(datalist)]
            elif fieldname == 'stp': datalist = [computestp(datalist)]
            # GSR in fields are T(K), mixing ratio (kg/kg), and surface pressure (Pa)
            elif fieldname == 'thetae': datalist = [compute_thetae(datalist)]
            elif fieldname == 'rh2m': datalist = [compute_rh(datalist)]
            elif fieldname == 'sspf': datalist = [ compute_sspf(datalist, file_list['upp']) ]
            #elif fieldname == 'pbmin': datalist = [ datalist[1] - datalist[0][:,0,:] ]
            elif fieldname == 'pbmin': datalist = [ datalist[1] - datalist[0] ] # CSS changed above line for GRIB2
            elif fieldname in ['thck1000-500', 'thck1000-850'] : datalist = [ datalist[1]*0.1 - datalist[0]*0.1 ] # CSS added for thicknesses

        datadict[f] = []
        for data in datalist:


          if is_precip_diff(fieldname):
            if debug:
                  print("Deriving accumulated precipitation. Subract ensemble at first time from ensemble at last time")
            for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
            # last and first time in the requested time range
            ensemble_at_last_time  = data[-ENS_SIZE:]
            ensemble_at_first_time = data[:ENS_SIZE]
            data = ensemble_at_last_time - ensemble_at_first_time

         # perform mean/max/variance/etc to reduce 3D array to 2D
          spatial_dimensions = data.shape[1:] # works for 1D meshes like MPAS and 2D grids like WRF
          ntimes = int(data.shape[0]/ENS_SIZE)
          if (fieldtype == 'mean'):  data = np.mean(data, axis=0)
          elif (fieldtype == 'pmm'): data = compute_pmm(data)
          elif (fieldtype == 'max'): data = np.amax(data, axis=0)
          elif (fieldtype == 'min'): data = np.amin(data, axis=0)
          elif (fieldtype == 'var'): data = np.std(data, axis=0)
          elif (fieldtype == 'maxstamp'):
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (ntimes,ENS_SIZE,*spatial_dimensions))
                data = np.nanmax(data, axis=0)
          elif (fieldtype == 'summean'):
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (ntimes,ENS_SIZE,*spatial_dimensions))
                data = np.nansum(data, axis=0)
                data = np.nanmean(data, axis=0)
          elif (fieldtype == 'maxmean'):
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (ntimes,ENS_SIZE,*spatial_dimensions))
                data = np.nanmax(data, axis=0)
                data = np.nanmean(data, axis=0)
          elif (fieldtype == 'summax'):
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (ntimes,ENS_SIZE,*spatial_dimensions))
                data = np.nansum(data, axis=0)
                data = np.nanmax(data, axis=0)
          elif (fieldtype[0:3] == 'mem'):
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (ntimes,ENS_SIZE,*spatial_dimensions))
                data = np.nanmax(data, axis=0)
                data = data[member-1,:] 
          elif (fieldtype in ['prob', 'neprob', 'problt', 'probgt', 'neprobgt', 'neproblt']):
                if fieldtype in ['prob', 'neprob', 'probgt', 'neprobgt']: data = (data>=thresh).astype('float')
                elif fieldtype in ['problt', 'neproblt']: data = (data<thresh).astype('float')
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (ntimes,ENS_SIZE,*spatial_dimensions))
                data = np.nanmax(data, axis=0)
                if fieldtype in ['neprob','neprobgt','neproblt']: 
                    junk = webPlot()
                    junk.domain = domain# define correct domain so .ibox can be correct
                    junk.loadMap() 
                    data2d = junk.latlonGrid(data[0,:])
                    data3d = np.tile(data2d, (ENS_SIZE,1,1))
                    for i in range(1,ENS_SIZE):
                        data2d = junk.latlonGrid(data[i,:])
                        data3d[i] = data2d
                    miles_to_km = 1.60934
                    roi = 25 * miles_to_km / junk.min_grid_spacing_km 
                    if debug: print("roi",roi)
                    data = compute_neprob(data3d, roi=int(roi), sigma=float(fields['sigma']), type='gaussian')
                else: data = np.nanmean(data, axis=0) 
                data = data+0.001 #hack to ensure that plot displays discrete prob values
          elif (fieldtype in ['prob3d']):
                data = (data>=thresh).astype('float')
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0)
                data = np.reshape(data, (ntimes,ENS_SIZE,*spatial_dimensions))
                data = compute_prob3d(data, roi=14, sigma=float(fields['sigma']), type='gaussian')
          if debug: print('field '+ fieldname+ ' has shape', data.shape, 'max', data.max(), 'min', data.min())

          print(data.max())
          #kernel = np.ones((7,7))
          #data = ndimage.filters.convolve(data, kernel/float(kernel.sum()))
          # attach data arrays for each type of field (e.g. { 'fill':[data], 'barb':[data,data] })
          datadict[f].append(data)

        fh.close()

    return (datadict, missing_list)


def is_precip_diff(s):
    # Is this string an accumulated precipitation difference field?
    # Assume it is if it starts with "precip-", ends with "hr"
    if s[0:7] == "precip-" and s[-2:] == "hr":
        return True
    return False


def readGrid(file_dir):
    f = Dataset(file_dir, 'r')
    lats   = f.variables['XLAT'][0,:]
    lons   = f.variables['XLONG'][0,:]
    f.close()
    return (lats,lons)

def readMPASVertices(ifile="/glade/p/mmm/parc/schwartz/MPAS/15-3km_mesh/init.nc"):
    # Used for regridding field from vertex to cell. latVertex and lonVertex not needed
    fh = Dataset(ifile, "r")
    nEdgesOnCell = fh.variables['nEdgesOnCell'][:]
    verticesOnCell = fh.variables['verticesOnCell'][:]
    fh.close()
    return (nEdgesOnCell, verticesOnCell)

def readGridMPAS():
    fh = Dataset("/glade/p/mmm/parc/schwartz/MPAS/15-3km_mesh/init.nc", "r")
    latCell = fh.variables['latCell'][:]
    lonCell = fh.variables['lonCell'][:]
    areaCell = fh.variables['areaCell'][:] # units m^2
    min_grid_spacing_km = 2. * np.sqrt(areaCell.min()/np.pi/1000/1000)
    # min_grid_spacing_km used for grid spacing of interpolated lat-lon grid.
    fh.close()
    latCell, lonCell = np.degrees(latCell), np.degrees(lonCell) #convert radians to degrees
    lonCell[lonCell >= 180] = lonCell[lonCell >= 180] - 360
    return (latCell, lonCell, min_grid_spacing_km)

def saveNewMap(domstr='CONUS', wrfout=None):
    # if domstr is not in the dictionary, then use provided wrfout to create new domain
    if domstr not in domains:
        fh = Dataset(wrfout, 'r')
        lats = fh.variables['XLAT'][0,:]
        lons = fh.variables['XLONG'][0,:]
        ll_lat, ll_lon, ur_lat, ur_lon = lats[0,0], lons[0,0], lats[-1,-1], lons[-1,-1]
        lat_1, lat_2, lon_0 = fh.TRUELAT1, fh.TRUELAT2, fh.STAND_LON 
        fig_width = 1080
        fh.close()
    # else assume domstr is in dictionary
    elif 'file' in domains[domstr]:
        fh = Dataset(domains[domstr]['file'], 'r')
        lats = fh.variables['XLAT'][0,:]
        lons = fh.variables['XLONG'][0,:]
        ll_lat, ll_lon, ur_lat, ur_lon = lats[0,0], lons[0,0], lats[-1,-1], lons[-1,-1]
        lat_1, lat_2, lon_0 = fh.TRUELAT1, fh.TRUELAT2, fh.STAND_LON 
        if 'fig_width' in domains[domstr]: fig_width = domains[domstr]['fig_width']
        else: fig_width = 1080
        fh.close()
    else:
        ll_lat, ll_lon, ur_lat, ur_lon = domains[domstr]['corners']
        fig_width = domains[domstr]['fig_width']
        lat_1, lat_2, lon_0 = 32.0, 46.0, -101.0
        if domstr=='NA':
            lon_0 = -115.0

    dpi = 90 
    fig = plt.figure(dpi=dpi)
    m = Basemap(projection='lcc', resolution='i', llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat, \
                lat_1=lat_1, lat_2=lat_2, lon_0=lon_0, area_thresh=1000)

    # compute height based on figure width, map aspect ratio, then add some vertical space for labels/colorbar
    fig_width  = fig_width/float(dpi)
    fig_height = fig_width*m.aspect + 0.93
    #fig_height = fig_width*m.aspect + 1.25
    figsize = (fig_width, fig_height)
    fig.set_size_inches(figsize)

    # place map 0.7" from bottom of figure, leave rest of 0.93" at top for title (needs to be in figure-relative coords)
    #x,y,w,h = 0.01, 0.8/float(fig_height), 0.98, 0.98*fig_width*m.aspect/float(fig_height) #too much padding at top
    x,y,w,h = 0.01, 0.7/float(fig_height), 0.98, 0.98*fig_width*m.aspect/float(fig_height)
    ax = fig.add_axes([x,y,w,h])
    for i in list(ax.spines.values()): i.set_linewidth(0.5)

    m.drawcoastlines(linewidth=0.5, ax=ax)
    m.drawstates(linewidth=0.25, ax=ax)
    m.drawcountries(ax=ax)
    # avoid this error
    # UnicodeDecodeError: 'utf-8' codec can't decode byte 0xf1 in position 2: invalid continuation byte
    #m.drawcounties(linewidth=0.1, color='gray', ax=ax)


    # load lat/lons
    lats, lons, min_grid_spacing_km = readGridMPAS()
    delta_deg = min_grid_spacing_km / 111
    if m.lonmin > 180 or m.lonmax > 180:
        lons[lons<0] = lons[lons<0] + 360.  # change -180-0 to 180-360 to match m.lonmin and m.lonmax
    # Replace m.lonmax with ur_lon? Otherwise, risk getting +179.999 when NW corner crosses the datetime
    nlon = int((m.lonmax - m.lonmin)/delta_deg)
    nlat = int((m.latmax - m.latmin)/delta_deg)
    if nlon > 1500:
        nlon = 1500
    if nlat > 1500:
        nlat = 1500
    lon2d, lat2d = np.meshgrid(np.linspace(m.lonmin, m.lonmax, nlon), np.linspace(m.latmin,m.latmax,nlat))
    # Convert to map coordinates instead of latlon to avoid the need to specify latlon=True in contour and barb methods.
    x2d, y2d = m(lon2d,lat2d)
    # ibox: subscripts within lat/lon box # only used to speed up 1-D array triangulation and plotting
    ibox = (m.lonmin-1 <= lons ) & (lons < m.lonmax+1) & (m.latmin-1 <= lats) & (lats < m.latmax+1)
    lons = lons[ibox]
    lats = lats[ibox]
    x, y = m(lons,lats)

    # use .filled() to avoid error about masked arrays
    vtx, wts = interp_weights(np.vstack((lons.filled(),lats.filled())).T,np.vstack((lon2d.flatten(), lat2d.flatten())).T)

    pickle.dump((fig,ax,m,lons,lats,min_grid_spacing_km,delta_deg,lon2d,lat2d,x2d,y2d,ibox,x,y,vtx,wts), open('/glade/work/ahijevyc/share/rt_ensemble/python_scripts/rt2015_%s.pk'%domstr, 'wb'))

def drawOverlay(domstr='CONUS'):
    ll_lat, ll_lon, ur_lat, ur_lon = domains[domstr]['corners']
    fig_width = domains[domstr]['fig_width']
    lat_1, lat_2, lon_0 = 32.0, 46.0, -101.0
    if domstr=='NA':
        lon_0 = -115.0
    dpi = 90

    fig = plt.figure(dpi=dpi)
    m = Basemap(projection='lcc', resolution='i', llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat, \
                lat_1=lat_1, lat_2=lat_2, lon_0=lon_0, area_thresh=1000)

    # compute height based on figure width, map aspect ratio, then add some vertical space for labels/colorbar
    fig_width  = fig_width/float(dpi)
    #fig_height = fig_width*m.aspect + 0.93
    fig_height = fig_width*m.aspect + 1.25
    figsize = (fig_width, fig_height)
    fig.set_size_inches(figsize)

    # place map 0.8" from bottom of figure, leave 0.45" at top for title (needs to be in figure-relative coords)
    x,y,w,h = 0.01, 0.8/float(fig_height), 0.98, 0.98*fig_width*m.aspect/float(fig_height)
    #x,y,w,h = 0.01, 0.7/float(fig_height), 0.98, 0.98*fig_width*m.aspect/float(fig_height)
    ax = fig.add_axes([x,y,w,h])

    #drawcounties doesnt work when called by itself? so have to drawcoastines first with lw=0
    m.drawcoastlines(linewidth=0, ax=ax)
    m.drawcounties(ax=ax)
    ax.axis('off')
    plt.savefig('overlay_counties_%s.png'%domstr, dpi=90, transparent=True)

def compute_pmm(ensemble):
    members = ensemble.shape[0]
    spatial_dimensions = ensemble.shape[1:]
    ens_mean = np.mean(ensemble, axis=0)
    ens_dist = np.sort(ensemble.flatten())[::-1]
    pmm = ens_dist[::members]

    ens_mean_index = np.argsort(ens_mean.flatten())[::-1]
    temp = np.empty_like(pmm)
    temp[ens_mean_index] = pmm

    temp = np.where(ens_mean.flatten() > 0, temp, 0.0)
    return temp.reshape(spatial_dimensions)

def compute_neprob(ensemble, roi=0, sigma=0.0, type='gaussian'):
    if len(ensemble.shape) < 3:
        print('compute_neprob: needs ensemble of 2D arrays, not 1D arrays')
        sys.exit(1)
    y,x = np.ogrid[-roi:roi+1, -roi:roi+1]
    kernel = x**2 + y**2 <= roi**2
    ens_roi = ndimage.filters.maximum_filter(ensemble, footprint=kernel[np.newaxis,:])

    ens_mean = np.nanmean(ens_roi, axis=0)
    #ens_mean = np.nanmean(ensemble, axis=0)

    if (type == 'uniform'):
        y,x = np.ogrid[-sigma:sigma+1, -sigma:sigma+1]
        kernel = x**2 + y**2 <= sigma**2
        ens_mean = ndimage.filters.convolve(ens_mean, kernel/float(kernel.sum()))
    elif (type == 'gaussian'):
        ens_mean = ndimage.filters.gaussian_filter(ens_mean, sigma)
    return ens_mean

def compute_prob3d(ensemble, roi=0, sigma=0.0, type='gaussian'):
    print(ensemble.shape)
    y,x = np.ogrid[-roi:roi+1, -roi:roi+1]
    kernel = x**2 + y**2 <= roi**2
    ens_roi = ndimage.filters.maximum_filter(ensemble, footprint=kernel[np.newaxis,np.newaxis,:])

    print(ens_roi.shape)
    ens_mean = np.nanmean(ens_roi, axis=1)
    print(ens_mean.shape)
    ens_mean = ndimage.filters.gaussian_filter(ens_mean, [2,20,20])
    return ens_mean[3,:]

def computestp(data):
    '''Compute STP with data array of [sbcape,sblcl,0-1srh,ushr06,vshr06]'''
    sbcape_term = (data[0]/1500.0)
    
    lcl_term = ((2000.0 - data[1])/1000.0)
    lcl_term = np.where(data[1] < 1000.0, 1.0, lcl_term)
    lcl_term = np.where(data[1] > 2000.0, 0.0, lcl_term)

    srh_term = (data[2]/150.0)

    ushear06 = data[3]-data[5]
    vshear06 = data[4]-data[6]
    shear06 = np.sqrt(ushear06**2 + vshear06**2) #this will be in knots (converted prior to fn)
    shear_term = (shear06/38.87)
    shear_term = np.where(shear06 > 58.32, 1.5, shear_term)
    shear_term = np.where(shear06 < 24.3, 0.0, shear_term)

    stp = (sbcape_term * lcl_term * srh_term * shear_term)
    # RS: this stopped working on 24 June 2016 - apparently stp not a masked array? replace with similar call, but may not be needed
    #stp = stp.filled(0.0) #fill missing values with 0s (apparently lcl_height missing along boundaries?)
    stp = np.ma.filled(stp, 0.0)
    return stp

def compute_sspf(data, cref_files):
    fh = MFDataset(cref_files)
    #cref = fh.variables['REFD_MAX'][:] #in diag files
    cref = fh.variables['REFL_MAX_COL'][:] #in upp files
    fh.close()

    # if all times are in one file, then need to reshape and extract desired times
    #cref = cref.reshape((10,49,cref.shape[1],cref.shape[2])) # reshape
    #cref = cref[:,13:36+1,:] # extract desired times
    #cref = np.swapaxes(cref,0,1) # flip first two axes so time is first
    #cref = cref.reshape((10*((13+1)-36)),cref.shape[2],cref.shape[3]) #reshape again

    # flag nearby points
    #kernel = np.ones((5,5))
    cref_mask = (cref>=30.0)
    #cref_mask = ndimage.filters.maximum_filter(cref_mask, footprint=kernel[np.newaxis,:])
 
    #wind_cref_hits = (data[1]>=25.0)
    wind_cref_hits = np.logical_and( (data[1]>=25.0), cref_mask)

    sspf = np.logical_or( np.logical_or((data[0]>=75), wind_cref_hits), (data[2]>=1))
    return sspf

def compute_uhratio(data):
    return np.where(data[1]>50, data[0]/data[1], 0.0)

def compute_thetae(data):
    # GSR constants for theta E calc
    P0 = 100000.0  # (Pa)
    Rd = 287.04    # (J/Kg)
    Cp = 1005.7    # (J/Kg)
    Lv = 2501000.0 # (J/Kg)
    LvoCp = Lv/Cp
    RdoCp = Rd/Cp
    return (((data[0]-32.0)/1.8)+273.15+LvoCp*data[1]) * ((P0/data[2])**RdoCp)

def compute_rh(data):
    t2 = (data[0]-32.0)/1.8 +  273.15 #temp in K
    psfc = data[1] #pres in Pa?
    q2 = data[2] #qvapor mixing ratio
    L_over_Rv = 5418.12

    es = 611.0 * np.exp(L_over_Rv*(1.0/273.0 - 1.0/t2))
    qsat = 0.622 * es / (psfc - es)
    rh = q2 / qsat
    return 100*rh

# from https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids
def interp_weights(xyz, uvw):
    tri = qhull.Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    d = xyz.shape[1]
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))


def showKeys():
    print(list(fieldinfo.keys()))
    sys.exit()

def computefrzdepth(t):
    frz_at_surface = np.where(t[0,:] < 33, True, False) #pts where surface T is below 33F
    max_column_t = np.amax(t, axis=0)
    above_frz_aloft = np.where(max_column_t > 32, True, False) #pts where max column T is above 32F

