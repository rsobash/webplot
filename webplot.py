import argparse
import cartopy
from collections import defaultdict
import datetime
from fieldinfo import domains, readNCLcm
import logging
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from metpy.units import units
import mpas
from mpas import fieldinfo, mesh_config, makeEnsembleList
import numpy as np
import pandas as pd
import pdb
import pickle
import os
import scipy.ndimage as ndimage
from scipy.interpolate import griddata
from scipy.spatial import qhull
import subprocess
import re
import sys
import time
# created with f2py3 -c mpas_vort_cell.f90 -m vert2cell
# Had to fix /glade/u/home/ahijevyc/miniconda3/envs/webplot/lib/python3.11/site-packages/numpy/f2py/src/fortranobject.c
# Moved int i declaration outside for loop. (line 707)
import vert2cell
import xarray

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

class webPlot:
    '''A class to plot data from NCAR ensemble'''
    def __init__(self):
        self.opts = parseargs()
        self.initdate = pd.to_datetime(self.opts['date'])
        self.ENS_SIZE = self.opts['ENS_SIZE']
        self.autolevels = self.opts['autolevels']
        self.debug = self.opts['debug']
        self.domain = self.opts['domain']
        self.fhr = self.opts['fhr']
        self.meshstr = self.opts['meshstr']
        self.nbarbs = self.opts['nbarbs']
        self.nlat_max = self.opts['nlat_max']
        self.nlon_max = self.opts['nlon_max']
        self.title = self.opts['title']

        self.get_mpas_mesh()
        self.loadMap()
        self.data, self.missing_members = readEnsemble(self)
        self.plotFields()
        self.plotTitleTimes()
        self.saveFigure()
 
    def createFilename(self):
        for f in ['fill', 'contour','barb']: # CSS added this for loop and everything in it
           prefx = self.meshstr
           if 'name' in self.opts[f]:
              if 'thresh' in self.opts[f]:
                 prefx = self.opts[f]['name']+'_'+self.opts[f]['ensprod']+'_'+str(self.opts[f]['thresh'])   # CSS
              else:
                 prefx = self.opts[f]['name']+'_'+self.opts[f]['ensprod'] # CSS
              break

        shr = min(self.fhr)
        ehr = max(self.fhr)
        if len(self.fhr) == 1:
            outfile = f"{prefx}_f{shr:03.0f}_{self.domain}.png"
        else: # CSS
            outfile = f"{prefx}_f{shr:03.0f}-f{ehr:03.0f}_{self.domain}.png"

        # create yyyymmddhh/domain/ directory if needed
        subdir_path = os.path.join(os.getenv('TMPDIR'), self.opts['date'], self.domain)
        if not os.path.isdir(subdir_path):
            logging.warning(f"webPlot.createFilename(): making new output directory {subdir_path}")
            os.makedirs(subdir_path)
        # prepend subdir_path to outfile.
        outfile = os.path.join(subdir_path, outfile)
        outfile = os.path.realpath(outfile)
        return outfile

    def loadMap(self):
        self.pk_file = os.path.join(os.getenv('TMPDIR'), f"{self.meshstr}_{self.domain}_{self.nlon_max}x{self.nlat_max}.pk")
        if not os.path.exists(self.pk_file):
            saveNewMap(self)
        logging.debug(f"loadMap: use old pickle file {self.pk_file}")
        (self.ax, self.extent, self.lons, self.lats, 
                self.lon2d, self.lat2d, self.x2d, self.y2d, self.ibox, self.x, self.y, self.vtx, 
                self.wts) = pickle.load(open(self.pk_file, 'rb'))

    def drawOverlay(self):
        domain=self.domain
        logging.info(f"drawOverlay: domain={domain}")
        self.ax.axis('off')
        self.ax.set_extent(self.extent, crs=self.ax.projection)
        plt.savefig(f'overlay_counties_{domain}.png', transparent=True)

    def get_mpas_mesh(self):
        path = self.opts["init_file"]
        mpas_mesh = xarray.open_dataset(path)
        lonCell = mpas_mesh['lonCell']
        lonCell = np.degrees(lonCell) #convert radians to degrees
        mpas_mesh["lonCell"] = lonCell
        mpas_mesh["latCell"] = np.degrees(mpas_mesh["latCell"]) #convert radians to degrees
        self.mpas_mesh = mpas_mesh

    def plotTitleTimes(self):
        fontdict = {'family':'monospace', 'size':12, 'weight':'bold'}

        # place title and times above corners of map
        x0, y1 = self.ax.transAxes.transform((0,1))
        x0, y0 = self.ax.transAxes.transform((0,0))
        x1, y1 = self.ax.transAxes.transform((1,1))
        self.ax.text(x0, y1+10, self.title, fontdict=fontdict, transform=None)

        fontdict = {'family':'monospace'}
        initstr, validstr = self.getInitValidStr()
        self.ax.text(1, 1, initstr+"\n"+validstr, fontdict=fontdict, horizontalalignment='right', 
                verticalalignment="bottom", transform=self.ax.transAxes)

        # Plot missing members 
        if len(self.missing_members):
            missing_members = sorted(set([ (x%self.ENS_SIZE)+1 for x in self.missing_members ])) #get member number from missing indices
            missing_members_string = ', '.join(str(x) for x in missing_members)
            self.ax.text(x1-5, y0+5, 'Missing member #s: %s'%missing_members_string, horizontalalignment='right')

    def plotFields(self):
        if 'fill' in self.data:
            if self.opts['fill']['ensprod'] == 'paintball': self.plotPaintball()
            elif self.opts['fill']['ensprod'].endswith("stamp"): self.plotStamp()
            else: self.plotFill()
        
        if 'contour' in self.data:
            if self.opts['contour']['ensprod'] == 'spaghetti': self.plotSpaghetti()
            elif self.opts['contour']['ensprod'].endswith('stamp'): self.plotStamp()
            else: self.plotContour()
        
        if 'barb' in self.data:
            assert not self.opts['barb']['ensprod'].endswith('stamp'), "TODO: postage stamp barbs"
            self.plotBarbs()
  
    def plotFill(self):
        if self.opts['fill']['name'] == 'crefuh': self.plotReflectivityUH(); return
        if self.autolevels:
            min, max = self.data['fill'][0].min(), self.data['fill'][0].max()
            levels = np.linspace(min, max, num=10)
            cmap = colors.ListedColormap(self.opts['fill']['colors'])
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
        if not self.opts['fill']['ensprod'].startswith('neprob'):
            logging.info(f"plotFill: latlonGrid")
            data = self.latlonGrid(data)
        if self.opts['fill']['name'] in ['avo500', 'vort500', 'pbmin']:
            logging.info(f"smooth {data.name}")
            # use .values to preserve DataArray attributes.
            data.values = ndimage.gaussian_filter(data, sigma=4)

        cs1 = self.ax.contourf(self.x2d, self.y2d, data, levels=levels, cmap=cmap, norm=norm, extend='max')
        label = f"{data.long_name} [{data.units}]"
        self.plotColorbar(cs1, levels, tick_labels, extend, extendfrac, label=label)

    def plotReflectivityUH(self):
        levels = self.opts['fill']['levels']
        cmap = colors.ListedColormap(self.opts['fill']['colors'])
        extend, extendfrac = 'neither', 0.0
        norm = colors.BoundaryNorm(levels, cmap.N)
        tick_labels = levels[:-1]

        logging.info(f"plotReflectivityUH: latlonGrid {self.data['fill'][0].long_name}")
        data = self.latlonGrid(self.data["fill"][0])
        cs1 = self.ax.contourf(self.x2d, self.y2d, data, levels=levels, cmap=cmap, norm=norm, extend='max')
        logging.info(f"plotReflectivityUH: latlonGrid {self.data['fill'][1].long_name}")
        uh =  self.latlonGrid(self.data['fill'][1])
        minUH = 100
        logging.info(f"plotReflectivityUH: UH shading and contour above {minUH}")
        self.ax.contourf(self.x2d, self.y2d, uh, levels=[minUH,1000], colors='black', alpha=0.3)
        cs2 = self.ax.contour( self.x2d, self.y2d, uh, levels=[minUH], colors='k', linewidths=0.5)
               
        label = f"{data.long_name} [{data.units}], shaded {uh.long_name} above ${minUH} {uh.units}$"
        self.plotColorbar(cs1, levels, tick_labels, extend, extendfrac, label=label)

    def plotColorbar(self, cs, levels, tick_labels, extend='neither', extendfrac=0.0, label=""):
        cb = plt.colorbar(cs, ax=self.ax, location="bottom", orientation='horizontal', extend=extend, shrink=0.65, 
                anchor=(1,1), panchor=(1,0), aspect=55, pad=0.03, extendfrac=extendfrac, ticks=tick_labels, label=label)
        cb.ax.xaxis.set_label_position('top')
        cb.outline.set_linewidth(0.5)

    def interpolatetri(self, values, vtx, wts):
        return np.einsum('nj,nj->n', np.take(values, vtx), wts)

    def latlonGrid(self, data):
        logging.debug(f"latlonGrid: {data.name}")
        data = data.metpy.dequantify() # Allow units to transfer to gridded array via attribute
        # apply ibox to data
        data = data.isel(nCells=self.ibox)
        if hasattr(self, "vtx") and hasattr(self, "wts"):
            logging.debug("latlonGrid: interpolatetri(vtx and wts)")
            # by using .values, Avoid interpolatetri ValueError: dimensions ('nCells',) must have the same length as the number of data dimensions, ndim=2
            data_gridded = self.interpolatetri(data.values, self.vtx, self.wts)
            data_gridded = np.reshape(data_gridded, self.lat2d.shape)
        else:
            logging.info("latlonGrid: interpolate with griddata()")
            data_gridded = griddata((self.lons, self.lats), data, (self.lon2d, self.lat2d), method='nearest')
        data_gridded = xarray.DataArray(data = data_gridded, coords = dict(lat=self.lat2d[:,0], lon=self.lon2d[0]), attrs=data.attrs)
        return data_gridded

    def plotContour(self):
        if self.opts['contour']['name'] in ['sbcinh','mlcinh']: linewidth, alpha = 0.5, 0.75
        else: linewidth, alpha = 1.5, 1.0
        data = self.data['contour'][0]

        if not self.opts['contour']['ensprod'].startswith('neprob'):
            logging.info(f"plotContour: latlonGrid")
            data = self.latlonGrid(data)

        if self.opts['contour']['name'] in ['t2-0c']: data.values = ndimage.gaussian_filter(data, sigma=2)
        else: data.values = ndimage.gaussian_filter(data, sigma=25)

        cs2 = self.ax.contour(self.x2d, self.y2d, data, levels=self.opts['contour']['levels'], colors='k', linewidths=linewidth, alpha=alpha)
        plt.clabel(cs2, fontsize='small', fmt='%i')

    def plotBarbs(self):

        skip = max([*self.x2d.shape, *self.y2d.shape])/self.nbarbs
        skip = int(skip)
        logging.debug(f"plotBarbs: nbarbs={self.nbarbs} skip={skip}")

        if self.opts['fill']['name'] == 'crefuh': alpha=0.5
        else: alpha=1.0

        logging.debug(f"plotBarbs: starting barbs {[x.name for x in self.data['barb']]}")
        # skip interval intended for 2-D fields
        u = self.latlonGrid(self.data['barb'][0])[::skip,::skip].values.flatten()
        v = self.latlonGrid(self.data['barb'][1])[::skip,::skip].values.flatten()
        x = self.lon2d[::skip,::skip].flatten()
        y = self.lat2d[::skip,::skip].flatten()
        # transform orients the barbs properly on projection
        logging.info(f"plotBarbs: barbs")
        cs2 = self.ax.barbs(x, y, u, v, alpha=alpha, length=5.5, linewidth=0.25, sizes={'emptybarb':0.05}, transform=cartopy.crs.PlateCarree())
    
    def plotPaintball(self):
       rects, labels = [], []
       colorlist = self.opts['fill']['colors']
       levels = self.opts['fill']['levels']
       for i in range(self.data['fill'][0].shape[0]):
           cs = self.ax.contourf(self.x, self.y, self.data['fill'][0][i,self.ibox], tri=True, levels=levels, colors=[colorlist[i%len(colorlist)]], alpha=0.5)
           rects.append(plt.Rectangle((0,0),1,1,fc=colorlist[i%len(colorlist)]))
           labels.append("member %d"%(i+1))

       plt.legend(rects, labels, ncol=5, loc='right', bbox_to_anchor=(1.0,-0.05), fontsize=11, \
                  frameon=False, borderpad=0.25, borderaxespad=0.25, handletextpad=0.2)

    def plotSpaghetti(self):
       proxy = []
       colorlist = self.opts['contour']['colors']
       levels = self.opts['contour']['levels']
       data = self.data['contour'][0]
       data.values = ndimage.gaussian_filter(data, sigma=[0,4,4])
       for i in range(self.data['contour'][0].shape[0]):
           cs = self.ax.contour(self.x, self.y, data[i,:], levels=levels, colors='k', alpha=0.6, linewidths=1, linestyles='solid')
           #proxy.append(plt.Rectangle((0,0),1,1,fc=colorlist[i]))
       #plt.legend(proxy, ["member %d"%i for i in range(1,11)], ncol=5, loc='right', bbox_to_anchor=(1.0,-0.05), fontsize=11, \
       #           frameon=False, borderpad=0.25, borderaxespad=0.25, handletextpad=0.2)

    def plotStamp(self):

        num_rows, num_columns = 3, 4
        fig, axs = plt.subplots(nrows=num_rows, ncols=num_columns,
                subplot_kw={'projection':self.ax.projection},
                figsize=(14,8))
        fig.subplots_adjust(bottom=0.01, top=0.99, left=0.01, right=0.99, wspace=0.01, hspace=0.01)

        levels = self.opts['fill']['levels']
        cmap = colors.ListedColormap(self.opts['fill']['colors'])
        norm = colors.BoundaryNorm(levels, cmap.N)

        for member, ax in enumerate(axs.flatten()):
            ax.add_feature(cartopy.feature.COASTLINE.with_scale('110m'), linewidth=0.25)
            ax.add_feature(cartopy.feature.BORDERS.with_scale('110m'), linewidth=0.25)
            ax.add_feature(cartopy.feature.STATES.with_scale('110m'), linewidth=0.05)
            ax.add_feature(cartopy.feature.LAKES.with_scale('110m'), edgecolor='k', linewidth=0.25, facecolor='k', alpha=0.05)
            ax.text(0.03,0.97,member+1,ha="left",va="top",bbox=dict(boxstyle="square",lw=0.5,fc="white"), transform=ax.transAxes)
            # plot, unless file that has fill field is missing, then skip
            if member not in self.missing_members and member < self.ENS_SIZE:
                data = self.latlonGrid(self.data['fill'][0].isel(mem=member))
                logging.debug(f"plotStamp: starting contourf with regridded array member {member}")
                cs1 = ax.contourf(self.x2d, self.y2d, data, levels=levels, cmap=cmap, norm=norm, extend='max')
                ax.set_extent(self.extent, crs=self.ax.projection)
            if member >= self.ENS_SIZE:
                fig.delaxes(ax)


        # use every other tick for large colortables, remove last tick label for both
        if self.opts['fill']['name'] in ['goesch3', 'goesch4', 't2', 'precipacc' ]: ticks = levels[:-1][::2] # CSS added precipacc       
        else: ticks = levels[:-1] 

        label = f"{data.long_name} [{data.units}]"
        # add colorbar to figure
        cax = fig.add_axes([0.51,0.3,0.48,0.02])
        cb = plt.colorbar(cs1, cax=cax, orientation='horizontal', ticks=ticks, extendfrac=0.0, label=label)
        cb.outline.set_linewidth(0.5)
        cb.ax.tick_params(labelsize=9)

        # add init/valid text
        fontdict = {'family':'monospace', 'size':13, 'weight':'bold'}
        initstr, validstr = self.getInitValidStr() 

        fig.text(0.51, 0.22, self.title, fontdict=fontdict, transform=fig.transFigure)
        pos = axs.flatten()[-2].get_position()
        fontdict = {'family':'monospace'}
        fig.text(0.51, pos.y0+0.3*pos.height, " "+initstr+"\n"+validstr, fontdict=fontdict, ha="left", 
                transform=fig.transFigure)

        # add NCAR logo and text below logo
        xo, yo = fig.transFigure.transform((0.51,pos.y0))
        fig.figimage(plt.imread('ncar.png'), xo=xo, yo=yo)
        #plt.text(x+10, y+5, 'ensemble.ucar.edu', fontdict={'size':9, 'color':'#505050'}, transform=None)


    def getInitValidStr(self):
       initstr  = self.initdate.strftime('Init: %a %Y-%m-%d %H UTC')
       shr = min(self.fhr)
       ehr = max(self.fhr)
       fmt = '%a %Y-%m-%d %H UTC'
       validstr = "Valid: " + (self.initdate+datetime.timedelta(hours=shr)).strftime(fmt)
       if ehr > shr: # range of valid times
            validstr += " - " + (self.initdate+datetime.timedelta(hours=ehr)).strftime(fmt)
       return initstr, validstr

    def saveFigure(self, transparent=False):
        outfile = self.createFilename()
        # place NCAR logo 57 pixels below bottom of map, then save image 
        if 'ensprod' in self.opts['fill']:  # CSS needed incase not a fill plot
            if not transparent and not self.opts['fill']['ensprod'].endswith('stamp'):
                img = plt.imread('ncar.png')
                self.ax.figure.figimage(img, xo=10, yo=10)
                #plt.text(x+10, y-54, 'ensemble.ucar.edu', fontdict={'size':9, 'color':'#505050'}, transform=None)

        self.ax.set_extent(self.extent, crs=self.ax.projection)

        plt.savefig(outfile, transparent=transparent, bbox_inches="tight")
        
        if self.opts['convert']:
            # reduce colors to shrink file size
            if not self.opts['fill']: ncolors = 48 #if no fill field exists
            elif "prob" in self.opts['fill']['ensprod']: ncolors = 48
            elif self.opts['fill']['name'] in ['crefuh']: ncolors = 48
            else: ncolors = 255
            if os.environ['NCAR_HOST'] == "cheyenne":
                quant = '/glade/u/home/ahijevyc/bin_cheyenne/pngquant'
            else:
                quant = '/glade/u/home/ahijevyc/bin/pngquant'
            beforesize = os.path.getsize(outfile)
            command = f"{quant} {ncolors} {outfile} --ext=.png --force"
            logging.info(f"{beforesize} bytes before reducing colors")
            logging.debug(command)
            ret = subprocess.check_call(command.split())
        plt.clf()
        fsize = os.path.getsize(outfile)
        logging.info(f"created {outfile} {fsize} bytes")

def parseargs():
    '''Parse arguments and return dictionary of fill, contour and barb field parameters'''

    parser = argparse.ArgumentParser(description='Web plotting script for NCAR ensemble')
    parser.add_argument('date', help='initialization datetime')
    parser.add_argument('--autolevels', action='store_true', help='use min/max to determine levels for plot')
    parser.add_argument('-b', '--barb', help='barb field (FIELD_PRODUCT_THRESH)')
    parser.add_argument('-c', '--contour', help='contour field (FIELD_PRODUCT_THRESH)')
    parser.add_argument('-con', '--convert', default=True, action='store_false', help='run final image through imagemagick')
    parser.add_argument('-d', '--debug', action='store_true', help='turn on debugging')
    parser.add_argument('--domain', type=str, choices=domains.keys(), default="CONUS", help='domain of plot')
    parser.add_argument('-f', '--fill', help='fill field (FIELD_PRODUCT_THRESH), FIELD options:'
            f"{','.join(list(fieldinfo.keys()))} PRODUCT may be one of [max,maxstamp,min,mean,meanstamp,"
            "prob,neprob,problt,neproblt,paintball,stamp,spaghetti]")
    parser.add_argument('--fhr', nargs='+', type=float, default=[12], help='list of forecast hours')
    parser.add_argument('--meshstr', type=str, default='uni', help='mesh id or path to defining mesh')
    parser.add_argument('--nbarbs', type=int, default=32, help='max barbs in one dimension')
    parser.add_argument('--nlon_max', type=int, default=1500, help='max pts in longitude dimension')
    parser.add_argument('--nlat_max', type=int, default=1500, help='max pts in latitude dimension')
    parser.add_argument('--sigma', default=2, help='smooth probabilities using gaussian smoother')
    parser.add_argument('-t', '--title', help='title for plot')

    opts = vars(parser.parse_args()) # argparse.Namespace in form of dictionary
    if opts["debug"]:
        logging.getLogger().setLevel(logging.DEBUG)

    # Based on meshstr, define ENS_SIZE, and init_file
    meshstr = opts["meshstr"]
    if meshstr in mesh_config:
        opts["ENS_SIZE"], opts["init_file"] = mesh_config[meshstr]
    else:
        assert os.path.exists(meshstr), (f"--meshstr must be a recognized mesh id {mesh_config.keys()} " 
                "or path to file (with lat/lon). not '{meshstr}'")
        opts["init_file"] = meshstr
        opts["ENS_SIZE"] = 1
        meshstr = os.path.basename(os.path.dirname(meshstr))

    # now, convert slash-delimited fill, contour, and barb args into dicts
    for f in ['contour','barb','fill']:
        thisdict = {}
        if opts[f] is not None:
            input = opts[f].split('/')

            assert len(input) > 1, f"{f} has 2-3 components separated by /. Add '/mean'?"

            thisdict['name']      = input[0]
            thisdict['ensprod']   = input[1]
            thisdict['arrayname'] = fieldinfo[input[0]]['fname'] # name of variable in input file
            
            # assign contour levels and colors
            if (input[1] in ['prob', 'neprob', 'probgt', 'problt', 'neprobgt', 'neproblt', 'prob3d']):
                if len(input)<3:
                    logging.error(f"your {f} option has less than 3 components. It needs name, ensprod, and thresh.")
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
                if input[0].startswith('hgt'):
                    thisdict['levels']  = [2,4,6,8,10,15,20,25,30,35,40,45,50,55,60,65,70,75] #hgt 
                    thisdict['colors']  = readNCLcm('wind_17lev')
                elif input[0].startswith('spe'):
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
            
        opts[f] = thisdict
    return opts

     
def readEnsemble(Plot):
    initdate = Plot.initdate
    fhr = Plot.fhr
    fields = Plot.opts
    ENS_SIZE = Plot.ENS_SIZE
    ''' Reads in desired fields and returns 2-D arrays of data for each field (barb/contour/field) '''
    logging.debug(fields)
    
    datadict = {}
    file_list = []
    missing_list = []
    if Plot.meshstr in ['15-3km_mesh', '15km_mesh']:
        file_list, missing_list = makeEnsembleList(Plot) #construct list of files
    elif Plot.meshstr == 'uni' and ENS_SIZE == 1:
        idir = "/glade/campaign/mmm/parc/ahijevyc/MPAS"
        file_list = [os.path.join(idir, Plot.meshstr, initdate.strftime('%Y%m%d%H'),
            (initdate+datetime.timedelta(hours=f)).strftime('diag.%Y-%m-%d_%H.%M.%S.nc')) for f in fhr] 
    else:
        logging.error("no ensemble files. Exiting.")
        logging.error("Perhaps add --mesh mpas to make_webplot.py command line")
        sys.exit(1)

    # loop through fill field, contour field, barb field and retrieve required data
    for f in ['fill', 'contour', 'barb']:
        if not list(fields[f].keys()): continue
        logging.info(f"read {fields[f]['name']}")
        
        # save some variables for use in this function
        arrays = fields[f]['arrayname']
        ensprod = fields[f]['ensprod']
        fieldname = fields[f]['name']
        if ensprod in ['prob', 'neprob', 'problt', 'probgt', 'neprobgt', 'neproblt', 'prob3d']: thresh = fields[f]['thresh']
        if ensprod.startswith('mem'): member = int(ensprod[3:])
        
        # open xarray Dataset
        logging.info(f"opening xarray Dataset. {len(fhr)} forecast hours x {ENS_SIZE} members")
        paths = file_list
        logging.debug(f"1-d paths = {paths}")
        paths = np.array(paths).reshape(len(fhr), ENS_SIZE).tolist()
        logging.debug(f"2-d paths = {paths}")
        # Use a nested list of paths so xarray assigns data to both "Time" and "mem" dimensions.
        fh = xarray.open_mfdataset(paths,engine='netcdf4',combine="nested", concat_dim=['Time','mem'])
        # turn XTIME = b'2017-05-02_00:00:00                                             ' to proper datetime.
        # use .load() or only a '2' will be read
        fh["Time"] = pd.to_datetime(fh.xtime.load().isel(mem=0).astype(str).str.strip(), format="%Y-%m-%d_%H:%M:%S") 
  
        # loop through each field, wind fields will have two fields that need to be read
        datalist = []
        for n,array in enumerate(arrays):
            logging.debug(f'Reading data {array}')
            #read in 3D array (times*members,ny,nx) from file object
            if 'arraylevel' in fields[f]:
                if isinstance(fields[f]['arraylevel'], list): level = fields[f]['arraylevel'][n]
                else: level = fields[f]['arraylevel']
            else: level = None
            
            if level == 'max':  data = fh[array].max(dim=level)
            elif level is None: data = fh[array]
            else:               data = fh[array].sel(level=level)

            if 'nVertices' in fh[array].dims:
                logging.info(f"{fieldname} is on vertices. Put on cells")
                # .load() to avoid dask PerformanceWarning: Reshaping is producing a large chunk.
                data = data.load().stack(TimeMem=("Time","mem")).T 
                mpas_mesh = Plot.mpas_mesh
                nEdgesOnCell = mpas_mesh.nEdgesOnCell
                verticesOnCell = mpas_mesh.verticesOnCell
                maxEdges = mpas_mesh.maxEdges.size
                nVertLevels = data.TimeMem.size # for now treat Time like vertical dimension
                nCells = mpas_mesh.nCells.size
                nVertices = data.nVertices.size
                # verticesOnCell is the transpose of what mpas_vort_cell1 expects
                dataCells = vert2cell.vert2cell(nEdgesOnCell, verticesOnCell.T, maxEdges, nVertLevels, nCells, nVertices, data)
                # Assign to new DataArray with nCells dimension. transfer DataArray attributes
                data = xarray.DataArray(data = dataCells, coords = dict(TimeMem=data.TimeMem, nCells=fh.nCells), attrs=data.attrs).unstack()

            # change units for certain fields
            if array in ['u10','v10']:                                  data = data.metpy.convert_units("knot")
            elif "zonal" in array or "meridional" in array:             data = data.metpy.convert_units("knot")
            elif array in ['dewpoint_surface', 't2m']:                  data = data.metpy.convert_units("degF")
            elif array in ['precipw','rainnc','grpl_max','SNOW_ACC_NC']:data = data.metpy.convert_units("inch")
            elif array in ['T_PL', 'TD_PL', 'SFC_LI']:                  data = data.metpy.convert_units("Celsius")
            elif array.startswith("temp"):                              data = data.metpy.convert_units("Celsius")
            elif array.startswith("dewp"):                              data = data.metpy.convert_units("Celsius")
            elif array == 'mslp' or 'PB' in array:                      data = data.metpy.convert_units("hPa")
            
            datalist.append(data)

        # these are derived fields, we don't have in any of the input files but we can compute
        if 'name' in fields[f]:
            if fieldname in ['shr06mag', 'shr01mag']:
                logging.info(f"derive {fieldname} from {arrays}")
                # Assume datalist is 4-element list: [bottom_u, bottom_v, top_u, top_v]
                shearmag = ((datalist[0]-datalist[2])**2 + (datalist[1]-datalist[3])**2)**0.5
                shearmag.attrs["long_name"] = fieldname.replace("shr","").replace("0","0-").replace("mag", "km shear magnitude")
                datalist = [shearmag]
            elif fieldname.startswith('speed') or fieldname == 'bunkmag':
                logging.info(f"derive speed from {arrays}")
                speed = (datalist[0]**2 + datalist[1]**2)**0.5
                speed.attrs["long_name"] = datalist[0].attrs["long_name"].replace("zonal wind"," wind")
                datalist = [speed]
            elif fieldname in ['shr06','shr01']: datalist = [datalist[0]-datalist[2], datalist[1]-datalist[3]]
            elif fieldname == 'uhratio': datalist = [compute_uhratio(datalist)]
            elif fieldname == 'stp': datalist = [computestp(datalist)]
            # GSR in fields are T(K), mixing ratio (kg/kg), and surface pressure (Pa)
            elif fieldname == 'thetae': datalist = [compute_thetae(datalist)]
            elif fieldname == 'rh2m': datalist = [compute_rh(datalist)]
            #elif fieldname == 'pbmin': datalist = [ datalist[1] - datalist[0][:,0,:] ]
            elif fieldname == 'pbmin': datalist = [ datalist[1] - datalist[0] ] # CSS changed above line for GRIB2
            elif fieldname in ['thck1000-500', 'thck1000-850'] : datalist = [ datalist[1]*0.1 - datalist[0]*0.1 ] # CSS added for thicknesses

        datadict[f] = []
        for data in datalist:
            if fieldname.startswith("precip"):
                fmt = '%a %Y-%m-%d %H UTC'
                logging.info(f"Derive accumulated precipitation {data.Time[0].dt.strftime(fmt).data} - {data.Time[-1].dt.strftime(fmt).data}")
                # subtract first time from last time
                ensemble_at_last_time  = data.isel(Time=[-1]) # brackets preserve Time dimension.
                ensemble_at_first_time = data.isel(Time=[0])
                data = ensemble_at_last_time # use last time for output data
                data.data -= ensemble_at_first_time.data # .data preserves quantity units, metadata, avoids pint.errors.DimensionalityError: Cannot convert from 'dimensionless' to 'inch'

            logging.info(f"get {ensprod} of {data.shape} data")
            if (ensprod == 'mean'):
                data = data.mean(dim=["Time","mem"], keep_attrs=True)
            elif (ensprod == 'pmm'): 
                data = compute_pmm(data)
            elif (ensprod == 'max'):
                data = data.max(dim=["Time","mem"], keep_attrs=True)
            elif (ensprod == 'min'):
                data = data.min(dim=["Time","mem"], keep_attrs=True)
            elif (ensprod == 'var'):
                data = data.std(dim=["Time","mem"], keep_attrs=True)
            elif (ensprod == 'maxstamp'):
                data = data.max(dim="Time", keep_attrs=True)
            elif (ensprod == 'meanstamp'):
                data = data.mean(dim="Time", keep_attrs=True)
            elif (ensprod == 'summean'):
                data = data.sum(dim="Time", keep_attrs=True)
                data = data.mean(dim="mem", keep_attrs=True)
            elif (ensprod == 'maxmean'):
                data = data.max(dim="Time", keep_attrs=True)
                data = data.mean(dim="mem", keep_attrs=True)
            elif (ensprod == 'summax'):
                data = data.sum(dim="Time", keep_attrs=True)
                data = data.max(dim="mem", keep_attrs=True)
            elif (ensprod[0:3] == 'mem'):
                data = data.sel(mem=member)
            elif 'prob' in ensprod:
                u = data.metpy.units
                if ensprod.endswith('lt'):
                    data.values = data.values < thresh
                    data.attrs["long_name"] = "less than {thresh} {u}"
                else:
                    data.values = data.values >= thresh
                    data.attrs["long_name"] = f"greater than or equal to {thresh} {u}"
                if ensprod.startswith("ne"):
                    # grid spacing of interpolated lat-lon grid.
                    grid_spacing = units.km * np.sqrt(Plot.mpas_mesh["areaCell"].values.min())/1000
                    nbrhd = 25 * units.miles
                    roi = nbrhd / grid_spacing
                    roi = roi.to_base_units()
                    logging.info(f"compute neighborhood probability with radius {roi:.2f}")
                    data = data.stack(TimeMem=("Time","mem")).groupby("TimeMem").apply(Plot.latlonGrid)
                    data = compute_neprob(data, roi=roi, sigma=float(Plot.opts['sigma']), type='gaussian')
                    data = data.unstack()
                    data.attrs["long_name"] += f" within {nbrhd:~}"
                data = data.mean(dim=["Time","mem"], keep_attrs=True)
                data.attrs["long_name"] = "probability " + data.attrs["long_name"]
                data.attrs["units"] = "dimensionless"
            elif (ensprod in ['prob3d']):
                data = (data.values>=thresh).astype('float')
                data = compute_prob3d(data, roi=14, sigma=float(Plot.opts['sigma']), type='gaussian')
            logging.debug(f'field {fieldname} has shape {data.shape} range {data.min()}-{data.min()}')

            datadict[f].append(data)

    return datadict, missing_list





def saveNewMap(plot):
    logging.info(f"saveNewMap: pk_file={plot.pk_file}")
    ll_lat, ll_lon, ur_lat, ur_lon = domains[plot.domain]['corners']
    lat_1, lat_2, lon_0 = 32., 46., -101.
    fig = plt.figure()
    
    proj = cartopy.crs.LambertConformal(central_longitude=lon_0, standard_parallels=(lat_1,lat_2))
    (llx, lly, llz), (urx, ury, urz) = proj.transform_points(
            cartopy.crs.PlateCarree(), 
            np.array([ll_lon, ur_lon]), 
            np.array([ll_lat, ur_lat])
            )
    ul_lon, ul_lat = cartopy.crs.PlateCarree().transform_point(llx, ury, proj)
    lr_lon, lr_lat = cartopy.crs.PlateCarree().transform_point(urx, lly, proj)
    lc_lon, lc_lat = cartopy.crs.PlateCarree().transform_point(0, lly, proj)
    uc_lon, uc_lat = cartopy.crs.PlateCarree().transform_point(0, ury, proj)

    # To save time, triangulate within a lat/lon bounding box.
    # Get extreme longitudes and latitudes within domain so the entire domain is covered.
    # These were attributes of Basemap object, but not cartopy.
    # Extreme longitudes are probably in upper left and upper right corners, but also check lower corners.
    lonmin = min([ll_lon, ul_lon])
    lonmax = max([lr_lon, ur_lon])
    # Extreme latitudes are probably in lower center and upper center (lc_lat, uc_lat).
    latmin = min([ll_lat, lc_lat, lr_lat])
    latmax = max([ul_lat, uc_lat, ur_lat])

    extent = (llx, urx, lly, ury) # in projection coordinates (x,y) meters

    # y/x aspect ratio was an attribute of Basemap object, but not cartopy.
    aspect = (ury - lly) / (urx - llx)

    # Constant figure width, no matter the aspect ratio.
    fig.set_size_inches(16,16*aspect)

    ax = plt.axes(projection = proj)
    for i in list(ax.spines.values()): i.set_linewidth(0.5)

    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'), linewidth=0.25)
    ax.add_feature(cartopy.feature.BORDERS.with_scale('50m'), linewidth=0.25)
    ax.add_feature(cartopy.feature.STATES.with_scale('50m'), linewidth=0.05)
    ax.add_feature(cartopy.feature.LAKES.with_scale('50m'), edgecolor='k', linewidth=0.25, facecolor='k', alpha=0.05)
    if plot.domain != "CONUS":
        logging.debug("draw counties")
        # Create custom cartopy feature COUNTIES that can be added to the axes.
        reader = cartopy.io.shapereader.Reader('/glade/work/ahijevyc/share/shapeFiles/cb_2013_us_county_500k/cb_2013_us_county_500k.shp')
        counties = list(reader.geometries())
        COUNTIES = cartopy.feature.ShapelyFeature(counties, cartopy.crs.PlateCarree())
        ax.add_feature(COUNTIES, facecolor="none", linewidth=0.1, alpha=0.25)

    # lat/lons from mpas_mesh file
    mpas_mesh = plot.mpas_mesh
    # min_grid_spacing is the grid spacing of interpolated lat-lon grid in meters.
    min_grid_spacing = np.sqrt(plot.mpas_mesh["areaCell"].values.min())

    lats = mpas_mesh["latCell"]
    lons = mpas_mesh["lonCell"]
    if lons.max() > 180:
        lons[lons>180] -= 360
    nlon = int((urx - llx)/min_grid_spacing) # calculate in projection space (meters)
    nlat = int((ury - lly)/min_grid_spacing)
    nlon = np.clip(nlon, 1, plot.nlon_max)
    nlat = np.clip(nlat, 1, plot.nlat_max)
    lon2d, lat2d = np.meshgrid(np.linspace(lonmin, lonmax, nlon), np.linspace(latmin,latmax,nlat))
    # Convert to map coordinates instead of latlon to avoid transform=PlateCarree in contour method.
    xyz = proj.transform_points(cartopy.crs.PlateCarree(), lon2d, lat2d)
    x2d = xyz[...,0]
    y2d = xyz[...,1]
    # ibox: subscripts within lat/lon box plus buffer. speed up triangulation in interp_weights
    ibox = (lonmin-1 <= lons ) & (lons < lonmax+1) & (latmin-1 <= lats) & (lats < latmax+1)
    lons = lons[ibox]
    lats = lats[ibox]
    xyz = proj.transform_points(cartopy.crs.PlateCarree(), lons, lats)
    x = xyz[...,0]
    y = xyz[...,1]

    logging.info(f"saveNewMap: triangulate {len(lons)} pts")
    vtx, wts = interp_weights(np.vstack((lons,lats)).T,np.vstack((lon2d.flatten(), lat2d.flatten())).T)

    pickle.dump((ax,extent,lons,lats,lon2d,lat2d,x2d,y2d,ibox,x,y,vtx,wts), open(plot.pk_file, 'wb'))
    
def compute_pmm(ensemble):
    members = ensemble.Time.size
    ens_mean = ensemble.mean(dim="Time", keep_attrs=True)
    logging.info(f"compute_pmm: sort {ensemble.values.size} ensemble values")
    ens_dist = np.sort(ensemble.values.flatten())[::-1]
    pmm = ens_dist[::members] # brilliant

    logging.debug(f"compute_pmm: sort {ens_mean.values.size} ens_mean values")
    ens_mean_index = np.argsort(ens_mean.values.flatten())[::-1]
    temp = np.empty_like(pmm)
    temp[ens_mean_index] = pmm

    temp = np.where(ens_mean.values.flatten() > 0, temp, 0.0)
    ens_mean.values = temp
    return ens_mean

def compute_neprob(ensemble, roi=0, sigma=0.0, type='gaussian'):
    roi = np.rint(roi) # round to nearest integer
    assert len(ensemble.dims) >= 3, 'compute_neprob: needs ensemble of 2D arrays, not 1D arrays'
    y,x = np.ogrid[-roi:roi+1, -roi:roi+1]
    kernel = x**2 + y**2 <= roi**2
    ens_roi = ndimage.filters.maximum_filter(ensemble, footprint=kernel[np.newaxis,:])

    if type == 'uniform':
        y,x = np.ogrid[-sigma:sigma+1, -sigma:sigma+1]
        kernel = x**2 + y**2 <= sigma**2
        smoothed = ndimage.filters.convolve(ens_roi, kernel/float(kernel.sum()))
    elif type == 'gaussian':
        # smooth last 2 dimensions (not TimeMember dimension)
        smoothed = ndimage.filters.gaussian_filter(ens_roi, (0,sigma,sigma))
    else:
        logging.error(f"compute_neprob: unknown filter {type}")
        sys.exit(1)
    ensemble.values = smoothed
    return ensemble

def compute_prob3d(ensemble, roi=0, sigma=0.):
    logging.info(f"compute_prob3d: roi={roi} sigma={sigma}")
    y,x = np.ogrid[-roi:roi+1, -roi:roi+1]
    kernel = x**2 + y**2 <= roi**2
    ens_roi = ndimage.filters.maximum_filter(ensemble, footprint=kernel[np.newaxis,np.newaxis,:])

    logging.info(f"ens_roi.shape={ens_roi.shape}")
    ens_mean = np.nanmean(ens_roi, axis=1)
    logging.info(f"ens_mean.shape={ens_mean.shape}")
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
    fh = xarray.open_mfdataset(cref_files)
    #cref = fh['REFD_MAX'] #in diag files
    cref = fh['REFL_MAX_COL'] #in upp files

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


def computefrzdepth(t):
    frz_at_surface = np.where(t[0,:] < 33, True, False) #pts where surface T is below 33F
    max_column_t = np.amax(t, axis=0)
    above_frz_aloft = np.where(max_column_t > 32, True, False) #pts where max column T is above 32F

if __name__ == "__main__":
    webPlot()

