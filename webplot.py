import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import *
from datetime import *
import pickle as pickle
import os, sys, time, argparse
import scipy.ndimage as ndimage
import subprocess
from fieldinfo import *
from netCDF4 import Dataset, MFDataset

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
        for f in ['fill', 'contour','barb']:
           if 'name' in self.opts[f]:
              if 'thresh' in self.opts[f]:
                 prefx = self.opts[f]['name']+'_'+self.opts[f]['ensprod']+'_'+str(self.opts[f]['thresh'])
              else:
                 prefx = self.opts[f]['name']+'_'+self.opts[f]['ensprod']
              break

        if self.shr == self.ehr:
            self.outfile = prefx+'_f'+'%03d'%self.shr+'_'+self.domain+'.png'
        else:
           self.outfile = prefx+'_f'+'%03d'%self.shr+'-f'+'%03d'%self.ehr+'_'+self.domain+'.png'

    def loadMap(self):
        # load pickle file containing figure and axes objects (should be pregenerated)
        PYTHON_SCRIPTS_DIR = os.getenv('PYTHON_SCRIPTS_DIR', '.')   
        self.fig, self.ax, self.m  = pickle.load(open('%s/%s.pk'%(PYTHON_SCRIPTS_DIR,self.domain), 'rb'))

        # get lat/lons from file here
        LATLON_FILE = os.getenv('LATLON_FILE', PYTHON_SCRIPTS_DIR+'/latlonfile.nc')
        self.lats, self.lons = readGrid(LATLON_FILE)
        self.x, self.y = self.m(self.lons,self.lats)

    def readEnsemble(self):
        self.data, self.missing_members = readEnsemble(self.initdate, timerange=[self.shr,self.ehr], fields=self.opts, debug=self.debug, ENS_SIZE=self.ENS_SIZE)

    def plotTitleTimes(self):
        fontdict = {'family':'monospace', 'size':12, 'weight':'bold'}

        # place title and times above corners of map
        x0, y1 = self.ax.transAxes.transform((0,1))
        x0, y0 = self.ax.transAxes.transform((0,0))
        x1, y1 = self.ax.transAxes.transform((1,1))
        self.ax.text(x0, y1+10, self.title, fontdict=fontdict, transform=None)

        initstr  = self.initdate.strftime('Init: %a %Y-%m-%d %H UTC') 
        if ((self.ehr - self.shr) == 0):
            validstr = (self.initdate+timedelta(hours=self.shr)).strftime('Valid: %a %Y-%m-%d %H UTC')
        else:
            validstr1 = (self.initdate+timedelta(hours=(self.shr-1))).strftime('%a %Y-%m-%d %H UTC')
            validstr2 = (self.initdate+timedelta(hours=self.ehr)).strftime('%a %Y-%m-%d %H UTC')
            validstr = "Valid: %s - %s"%(validstr1, validstr2)

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
       
    def plotFill(self):
        if self.opts['fill']['name'] == 'ptype': self.plotFill_ptype(); return
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

        # smooth some of the fill fields
        if self.opts['fill']['name'] == 'avo500': self.data['fill'][0] = ndimage.gaussian_filter(self.data['fill'][0], sigma=4)
        if self.opts['fill']['name'] == 'pbmin': self.data['fill'][0] = ndimage.gaussian_filter(self.data['fill'][0], sigma=2)
        
        cs1 = self.m.contourf(self.x, self.y, self.data['fill'][0], levels=levels, cmap=cmap, norm=norm, extend='max', ax=self.ax)

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

    def plotReflectivityUH(self):
        levels = self.opts['fill']['levels']
        cmap = colors.ListedColormap(self.opts['fill']['colors'])
        norm = colors.BoundaryNorm(levels, cmap.N)
        tick_labels = levels[:-1]

        cs1 = self.m.contourf(self.x, self.y, self.data['fill'][0], levels=levels, cmap=cmap, norm=norm, extend='max', ax=self.ax)
        self.m.contourf(self.x, self.y, self.data['fill'][1], levels=[75,1000], colors='black', ax=self.ax, alpha=0.3)
        self.m.contour(self.x, self.y, self.data['fill'][1], levels=[75], colors='k', linewidth=0.5, ax=self.ax)
               
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

    def plotContour(self):
        if self.opts['contour']['name'] in ['t2-0c']: data = ndimage.gaussian_filter(self.data['contour'][0], sigma=2)
        else: data = ndimage.gaussian_filter(self.data['contour'][0], sigma=10)

        if self.opts['contour']['name'] in ['sbcinh','mlcinh']: linewidth, alpha = 0.5, 0.75
        else: linewidth, alpha = 1.5, 1.0
        cs2 = self.m.contour(self.x, self.y, data, levels=self.opts['contour']['levels'], colors='k', linewidths=linewidth, ax=self.ax, alpha=alpha)
        plt.clabel(cs2, fontsize='small', fmt='%i')

    def plotBarbs(self):
        skip = self.opts['barb']['skip']
        if self.domain != 'CONUS': skip = 20

        if self.opts['fill']['name'] == 'crefuh': alpha=0.5
        else: alpha=1.0

        cs2 = self.m.barbs(self.x[::skip,::skip], self.y[::skip,::skip], self.data['barb'][0][::skip,::skip], self.data['barb'][1][::skip,::skip], \
                     color='black', alpha=alpha, length=5.5, linewidth=0.25, sizes={'emptybarb':0.05}, ax=self.ax)
    
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
           cs = self.m.contourf(self.x, self.y, self.data['fill'][0][i,:], levels=levels, colors=[colorlist[i%len(colorlist)]], ax=self.ax, alpha=0.5)
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

    def plotStamp(self):
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

               #print 'member', member, 'creating axes at', x, y
               thisax = fig.add_axes([x,y,w,h])

               thisax.axis('on')
               for axis in ['top','bottom','left','right']: thisax.spines[axis].set_linewidth(0.5)
               self.m.drawcoastlines(ax=thisax, linewidth=0.3)
               self.m.drawstates(linewidth=0.15, ax=thisax)
               self.m.drawcountries(ax=thisax, linewidth=0.3)
               thisax.text(0.03,0.97,member+1,ha="left",va="top",bbox=dict(boxstyle="square",lw=0.5,fc="white"), transform=thisax.transAxes)
               
               # plot, unless file that has fill field is missing, then skip
               if member not in self.missing_members[filename] and member < self.ENS_SIZE:
                   cs1 = self.m.contourf(self.x, self.y, self.data['fill'][0][memberidx,:], levels=levels, cmap=cmap, norm=norm, extend='max', ax=thisax)
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
       initstr  = self.initdate.strftime(' Init: %a %Y-%m-%d %H UTC')
       if ((self.ehr - self.shr) == 0):
            validstr = (self.initdate+timedelta(hours=self.shr)).strftime('Valid: %a %Y-%m-%d %H UTC')
       else:
            validstr1 = (self.initdate+timedelta(hours=(self.shr-1))).strftime('%a %Y-%m-%d %H UTC')
            validstr2 = (self.initdate+timedelta(hours=self.ehr)).strftime('%a %Y-%m-%d %H UTC')
            validstr = "Valid: %s - %s"%(validstr1, validstr2)

       fig.text(0.51, 0.22, self.title, fontdict=fontdict, transform=fig.transFigure)
       fig.text(0.51, 0.22 - 25/float(fig_height_px), initstr, transform=fig.transFigure)
       fig.text(0.51, 0.22 - 40/float(fig_height_px), validstr, transform=fig.transFigure)

       # add logo and text below logo
       x, y = fig.transFigure.transform((0.51,0))
       fig.figimage(plt.imread('ncar.png'), xo=x, yo=y+15, zorder=1000)
       plt.text(x+10, y+5, 'ensemble.ucar.edu', fontdict={'size':9, 'color':'#505050'}, transform=None)

    def saveFigure(self, trans=False):
        # place NCAR logo 57 pixels below bottom of map, then save image 
        if 'ensprod' in self.opts['fill']:  # CSS needed incase not a fill plot
           if not trans and self.opts['fill']['ensprod'] not in ['stamp', 'maxstamp']:
             x, y = self.ax.transAxes.transform((0,0))
             #self.fig.figimage(plt.imread('ncar.png'), xo=x, yo=(y-44), zorder=1000)
             #plt.text(x+10, y-54, 'ensemble.ucar.edu', fontdict={'size':9, 'color':'#505050'}, transform=None)

        plt.savefig(self.outfile, dpi=90, transparent=trans)
        
        if self.opts['convert']:
            #command = 'convert -colors 255 %s %s'%(self.outfile, self.outfile)
            if not self.opts['fill']: ncolors = 48 #if no fill field exists
            elif self.opts['fill']['ensprod'] in ['prob', 'neprob', 'probgt', 'problt', 'neprobgt', 'neproblt']: ncolors = 48
            elif self.opts['fill']['name'] in ['crefuh']: ncolors = 48
            else: ncolors = 255
            command = '/glade/u/home/sobash/pngquant/pngquant %d %s --ext=.png --force'%(ncolors,self.outfile)
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
    parser.add_argument('-sig', '--sigma', default=2, help='smooth probabilities using gaussian smoother')
    parser.add_argument('--debug', action='store_true', help='turn on debugging')

    opts = vars(parser.parse_args())
    
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
                thisdict['thresh']  = float(input[2])
                if int(opts['sigma']) != 40: thisdict['levels']  = np.arange(0.1,1.1,0.1)
                else: thisdict['levels']  = [0.02,0.05,0.1,0.15,0.2,0.25,0.35,0.45,0.6]
                thisdict['colors']  = readNCLcm('perc2_9lev')
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
                upp    = '%s/%s/post_rundir/mem_%d/fhr_%d/WRFTWO%02d.nc'%(EXP_DIR,yyyymmddhh,mem,hr,hr)
                if os.path.exists(wrfout): file_list['wrfout'].append(wrfout)
                else: missing_list['wrfout'].append(missing_index)
                if os.path.exists(diag): file_list['diag'].append(diag)
                else: missing_list['diag'].append(missing_index)
                if os.path.exists(upp): file_list['upp'].append(upp)
                else: missing_list['upp'].append(missing_index)
                missing_index += 1
    return (file_list, missing_list)

def readEnsemble(wrfinit, timerange=None, fields=None, debug=False, ENS_SIZE=10):
    ''' Reads in desired fields and returns 2-D arrays of data for each field (barb/contour/field) '''
    if debug: print(fields)
    
    datadict = {}
    file_list, missing_list = makeEnsembleList(wrfinit, timerange, ENS_SIZE) #construct list of files
 
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
        if debug: print(file_list[filename]) 
        fh = MFDataset(file_list[filename])
       
        # loop through each field, wind fields will have two fields that need to be read
        datalist = []
        for n,array in enumerate(arrays):
            if debug: print('Reading', array)

            #read in 3D array (times*members,ny,nx) from file object
            if 'arraylevel' in fields[f]:
                if isinstance(fields[f]['arraylevel'], list): level = fields[f]['arraylevel'][n]
                else: level = fields[f]['arraylevel']
            else: level = None
            
            if level == 'max':  data = np.amax(fh.variables[array][:,:,:,:], axis=1)
            elif level is None: data = fh.variables[array][:,:,:]
            else:               data = fh.variables[array][:,level,:,:]

            # change units for certain fields
            if array in ['U_PL', 'V_PL', 'UBSHR6','VBSHR6','UBSHR1', 'VBSHR1', 'U10','V10', 'U_COMP_STM', 'V_COMP_STM','S_PL','U_COMP_STM_6KM','V_COMP_STM_6KM']:  data = data*1.93 # m/s > kt
            elif array in ['DEWPOINT_2M', 'T2', 'AFWA_WCHILL', 'AFWA_HEATIDX']:   data = (data - 273.15)*1.8 + 32.0 # K > F 
            elif array in ['PREC_ACC_NC', 'PREC_ACC_C', 'AFWA_PWAT', 'PWAT', 'AFWA_RAIN', 'AFWA_SNOWFALL', 'AFWA_SNOW', 'AFWA_ICE', 'AFWA_FZRA','AFWA_RAIN_HRLY','AFWA_ICE_HRLY','AFWA_SNOWFALL_HRLY', 'AFWA_FZRA_HRLY']:   data = data*0.0393701 # mm > in             
            elif array in ['RAINNC', 'GRPL_MAX', 'SNOW_ACC_NC', 'AFWA_HAIL', 'HAILCAST_DIAM_MAX']:   data = data*0.0393701 # mm > in 
            elif array in ['T_PL', 'TD_PL', 'SFC_LI']:             data = data - 273.15 # K > C
            elif array in ['AFWA_MSLP', 'MSLP']:                   data = data*0.01 # Pa > hPa
            elif array in ['ECHOTOP']:                             data = data*3.28084# m > ft
            elif array in ['UP_HELI_MIN']:                         data = np.abs(data)
            elif array in ['AFWA_VIS', 'VISIBILITY']:              data = (data*0.001)/1.61  # m > mi
            elif array in ['SBCINH', 'MLCINH', 'W_DN_MAX']:        data = data*-1.0 # make cin positive
            elif array in ['PVORT_320K']:                          data = data*1000000 # multiply by 1e6
            elif array in ['SBT123_GDS3_NTAT','SBT124_GDS3_NTAT','GOESE_WV','GOESE_IR']: data = data -273.15 # K -> C
            elif array in ['HAIL_MAXK1', 'HAIL_MAX2D']:            data = data*39.3701 #  m -> inches
            elif array in ['PBMIN', 'PBMIN_SFC', 'BESTPBMIN', 'MLPBMIN', 'MUPBMIN']:     data = data*0.01 #  Pa -> hPa
#            elif array in ['LTG1_MAX1', 'LTG2_MAX', 'LTG3_MAX']:   data = data*0.20 #  scale down excess values
            
            datalist.append(data)

        # these are derived fields, we don't have in any of the input files but we can compute
        print(datalist[0].shape)
        if 'name' in fields[f]:
            if fieldname in ['shr06mag', 'shr01mag', 'bunkmag','speed10m']: datalist = [np.sqrt(datalist[0]**2 + datalist[1]**2)]
            elif fieldname == 'stp': datalist = [computestp(datalist)]
            # GSR in fields are T(K), mixing ratio (kg/kg), and surface pressure (Pa)
            elif fieldname == 'thetae': datalist = [compute_thetae(datalist)]
            elif fieldname == 'rh2m': datalist = [compute_rh(datalist)]
            #elif fieldname == 'pbmin': datalist = [ datalist[1] - datalist[0][:,0,:] ]
            elif fieldname == 'pbmin': datalist = [ datalist[1] - datalist[0] ] # CSS changed above line for GRIB2
            elif fieldname in ['thck1000-500', 'thck1000-850'] : datalist = [ datalist[1]*0.1 - datalist[0]*0.1 ] # CSS added for thicknesses
            elif fieldname == 'winter': datalist = [datalist[1] + datalist[2] + datalist[3]]

        datadict[f] = []
        for data in datalist:
          # perform mean/max/variance/etc to reduce 3D array to 2D
          if (fieldtype == 'mean'):  data = np.mean(data, axis=0)
          elif (fieldtype == 'pmm'): data = compute_pmm(data)
          elif (fieldtype == 'max'): data = np.amax(data, axis=0)
          elif (fieldtype == 'min'): data = np.amin(data, axis=0)
          elif (fieldtype == 'var'): data = np.std(data, axis=0)
          elif (fieldtype == 'maxstamp'):
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (data.shape[0]/ENS_SIZE,ENS_SIZE,data.shape[1],data.shape[2]))
                data = np.nanmax(data, axis=0)
          elif (fieldtype == 'summean'):
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (data.shape[0]/ENS_SIZE,ENS_SIZE,data.shape[1],data.shape[2]))
                data = np.nansum(data, axis=0)
                data = np.nanmean(data, axis=0)
          elif (fieldtype == 'maxmean'):
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (data.shape[0]/ENS_SIZE,ENS_SIZE,data.shape[1],data.shape[2]))
                data = np.nanmax(data, axis=0)
                data = np.nanmean(data, axis=0)
          elif (fieldtype == 'summax'):
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (data.shape[0]/ENS_SIZE,ENS_SIZE,data.shape[1],data.shape[2]))
                data = np.nansum(data, axis=0)
                data = np.nanmax(data, axis=0)
          elif (fieldtype[0:3] == 'mem'):
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (data.shape[0]/ENS_SIZE,ENS_SIZE,data.shape[1],data.shape[2]))
                print(fieldname)
                if fieldname in ['precip', 'precipacc']:
                    print('where we should be')
                    data = np.nanmax(data, axis=0)
                else:                      data = np.nanmax(data, axis=0)
                data = data[member-1,:] 
          elif (fieldtype in ['prob', 'neprob', 'problt', 'probgt', 'neprobgt', 'neproblt']):
                if fieldtype in ['prob', 'neprob', 'probgt', 'neprobgt']: data = (data>=thresh).astype('float')
                elif fieldtype in ['problt', 'neproblt']: data = (data<thresh).astype('float')
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0) #insert nan for missing files
                data = np.reshape(data, (data.shape[0]/ENS_SIZE,ENS_SIZE,data.shape[1],data.shape[2]))
                data = np.nanmax(data, axis=0)
                if (fieldtype in ['neprob', 'neproblt', 'neprobgt']): data = compute_neprob(data, roi=14, sigma=float(fields['sigma']), type='gaussian')
                else: data = np.nanmean(data, axis=0) 
                data = data+0.001 #hack to ensure that plot displays discrete prob values
          elif (fieldtype in ['prob3d']):
                data = (data>=thresh).astype('float')
                for i in missing_list[filename]: data = np.insert(data, i, np.nan, axis=0)
                data = np.reshape(data, (data.shape[0]/ENS_SIZE,ENS_SIZE,data.shape[1],data.shape[2]))
                data = compute_prob3d(data, roi=14, sigma=float(fields['sigma']), type='gaussian')
          if debug: print('field', fieldname, 'has shape', data.shape, 'max', data.max(), 'min', data.min())

          # attach data arrays for each type of field (e.g. { 'fill':[data], 'barb':[data,data] })
          datadict[f].append(data)

        fh.close()

    return (datadict, missing_list)

def readGrid(file_dir):
    f = Dataset(file_dir, 'r')
    lats   = f.variables['XLAT'][0,:]
    lons   = f.variables['XLONG'][0,:]
    f.close()
    return (lats,lons)

def saveNewMap(domstr='CONUS', wrfout=None):
    # THIS GENERATES A NEW PICKLE FILE FOR A DOMAIN (need to run once for a new domain)
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
    for i in ax.spines.values(): i.set_linewidth(0.5)

    m.drawcoastlines(linewidth=0.5, ax=ax)
    m.drawstates(linewidth=0.25, ax=ax)
    m.drawcountries(ax=ax)
    m.drawcounties(linewidth=0.1, ax=ax)

    pickle.dump((fig,ax,m), open('rt2015_%s.pk'%domstr, 'wb'))

def compute_pmm(ensemble):
    mem, dy, dx = ensemble.shape
    ens_mean = np.mean(ensemble, axis=0)
    ens_dist = np.sort(ensemble.flatten())[::-1]
    pmm = ens_dist[::mem]

    ens_mean_index = np.argsort(ens_mean.flatten())[::-1]
    temp = np.empty_like(pmm)
    temp[ens_mean_index] = pmm

    temp = np.where(ens_mean.flatten() > 0, temp, 0.0)
    return temp.reshape((dy,dx))

def compute_neprob(ensemble, roi=0, sigma=0.0, type='gaussian'):
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

    shear06 = np.sqrt(data[3]**2 + data[4]**2) #this will be in knots (converted prior to fn)
    shear_term = (shear06/38.87)
    shear_term = np.where(shear06 > 58.32, 1.5, shear_term)
    shear_term = np.where(shear06 < 24.3, 0.0, shear_term)

    stp = (sbcape_term * lcl_term * srh_term * shear_term)
    # RS: this stopped working on 24 June 2016 - apparently stp not a masked array? replace with similar call, but may not be needed
    #stp = stp.filled(0.0) #fill missing values with 0s (apparently lcl_height missing along boundaries?)
    stp = np.ma.filled(stp, 0.0)
    return stp

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
