import argparse
from collections import defaultdict
import datetime
from fieldinfo import domains
import logging
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import mpas
from mpas import fieldinfo, mesh_config, makeEnsembleList
import vert2cell # created with f2py3 -c mpas_vort_cell.f90 -m vert2cell
from mpl_toolkits.basemap import *
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
import xarray

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

class webPlot:
    '''A class to plot data from NCAR ensemble'''
    def __init__(self, domain=None):
        self.opts = parseargs()
        self.initdate = pd.to_datetime(self.opts['date'])
        self.ENS_SIZE = self.opts['ENS_SIZE']
        self.autolevels = self.opts['autolevels']
        self.debug = self.opts['debug']
        self.domain = domain
        self.fhr = self.opts['fhr']
        self.meshstr = self.opts['meshstr']
        self.nbarbs = self.opts['nbarbs']
        self.nlat_max = self.opts['nlat_max']
        self.nlon_max = self.opts['nlon_max']
        self.title = self.opts['title']

        self.get_mpas_mesh()
        self.createFilename()
        self.loadMap()
        self.data, self.missing_members = readEnsemble(self)
        self.plotFields()
        self.plotTitleTimes()
 
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
            self.outfile = f"{prefx}_f{shr:03.0f}_{self.domain}.png"
        else: # CSS
            self.outfile = f"{prefx}_f{shr:03.0f}-f{ehr:03.0f}_{self.domain}.png"

        # create yyyymmddhh/domain/ directory if needed
        subdir_path = os.path.join(os.getenv('TMPDIR'), self.opts['date'], self.domain)
        if not os.path.isdir(subdir_path):
            logging.warning(f"webPlot.createFilename(): making new output directory {subdir_path}")
            os.makedirs(subdir_path)
        # prepend subdir_path to outfile.
        self.outfile = os.path.join(subdir_path, self.outfile)
        self.outfile = os.path.realpath(self.outfile)

    def loadMap(self, overlay=False):
        pk_file = os.path.join(os.getenv('TMPDIR'), f"{self.meshstr}_{self.domain}_{self.nlon_max}x{self.nlat_max}.pk")
        if not os.path.exists(pk_file):
            saveNewMap(self, pk_file)
        logging.debug(f"loadMap: use old pickle file {pk_file}")
        (self.fig, self.ax, self.m, self.lons, self.lats, self.delta_deg,
                self.lon2d, self.lat2d, self.x2d, self.y2d, self.ibox, self.x, self.y, self.vtx, 
                self.wts) = pickle.load(open(pk_file, 'rb'))

    def get_mpas_mesh(self):
        path = self.opts["init_file"]
        mpas_mesh = xarray.open_dataset(path)
        lonCell = mpas_mesh['lonCell']
        lonCell = np.degrees(lonCell) #convert radians to degrees
        mpas_mesh["lonCell"] = lonCell
        mpas_mesh["latCell"] = np.degrees(mpas_mesh["latCell"]) #convert radians to degrees
        self.mpas_mesh = mpas_mesh

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
            assert not self.opts['contour']['ensprod'].endswith('stamp'), "TODO: postage stamp barbs"
            self.plotBarbs()
  
    def plotFill(self):
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
        data = self.latlonGrid(data)
        # regrid 1D mesh that needs to be smoothed
        if self.opts['fill']['name'] in ['avo500', 'vort500', 'pbmin']:
            # smooth some of the fill fields. use .values to preserve DataArray attributes.
            data.values = ndimage.gaussian_filter(data, sigma=4)
        elif self.opts['fill']['ensprod'].startswith('neprob'):
            miles_to_km = 1.60934
            roi = 25 * miles_to_km / self.min_grid_spacing_km 
            logging.debug(f"roi {roi}")
            data = compute_neprob(data, roi=int(roi), sigma=float(fields['sigma']), type='gaussian')
            data = data+0.001 #hack to ensure that plot displays discrete prob values

        cs1 = self.m.contourf(self.x2d, self.y2d, data, levels=levels, cmap=cmap, norm=norm, extend='max', ax=self.ax)
        label = f"{data.long_name} [{data.units}]"
        self.plotColorbar(cs1, levels, tick_labels, extend, extendfrac, label=label)

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
               
        self.plotColorbar(cs1, levels, tick_labels)

    def plotColorbar(self, cs, levels, tick_labels, extend='neither', extendfrac=0.0, label=""):
        # make axes for colorbar, 175px to left and 30px down from bottom of map 
        x0, y0 = self.ax.transAxes.transform((0,0))
        x, y = self.fig.transFigure.inverted().transform((x0+175,y0-29.5))
        cax = self.fig.add_axes([x,y,0.985-x,y/3.0])
        cb = plt.colorbar(cs, cax=cax, orientation='horizontal', extend=extend, extendfrac=extendfrac, ticks=tick_labels, label=label)
        cb.ax.xaxis.set_label_position('top')
        cb.outline.set_linewidth(0.5)

    def interpolatetri(self, values, vtx, wts):
        return np.einsum('nj,nj->n', np.take(values, vtx), wts)

    def latlonGrid(self, data):
        data = data.metpy.dequantify() # Allow units to transfer to gridded array via attribute
        # apply ibox to data
        data = data[self.ibox]
        if hasattr(self, "vtx") and hasattr(self, "wts"):
            logging.debug("latlonGrid: interpolatetri(vtx and wts)")
            # by using .values, Avoid interpolatetri ValueError: dimensions ('nCells',) must have the same length as the number of data dimensions, ndim=2
            data_gridded = self.interpolatetri(data.values, self.vtx, self.wts)
            data_gridded = np.reshape(data_gridded, self.lat2d.shape)
        else:
            logging.info("latlonGrid: interpolate to latlon grid with griddata()")
            data_gridded = griddata((self.lons, self.lats), data, (self.lon2d, self.lat2d), method='nearest')
        data_gridded = xarray.DataArray(data = data_gridded, coords = dict(lat=self.lat2d[:,0], lon=self.lon2d[0]), attrs=data.attrs)
        return data_gridded

    def plotContour(self):
        if self.opts['contour']['name'] in ['sbcinh','mlcinh']: linewidth, alpha = 0.5, 0.75
        else: linewidth, alpha = 1.5, 1.0
        data = self.data['contour'][0]
        data_gridded = self.latlonGrid(data)

        if self.opts['contour']['name'] in ['t2-0c']: data_gridded.values = ndimage.gaussian_filter(data_gridded, sigma=2)
        else: data_gridded.values = ndimage.gaussian_filter(data_gridded, sigma=25)

        cs2 = self.m.contour(self.x2d, self.y2d, data_gridded, levels=self.opts['contour']['levels'], colors='k', linewidths=linewidth, ax=self.ax, alpha=alpha)
        plt.clabel(cs2, fontsize='small', fmt='%i')

    def plotBarbs(self):

        skip = max([*self.x2d.shape, *self.y2d.shape])/self.nbarbs
        skip = int(skip)
        logging.debug(f"plotBarbs: nbarbs={self.nbarbs} skip={skip}")

        if self.opts['fill']['name'] == 'crefuh': alpha=0.5
        else: alpha=1.0

        logging.debug(f"plotBarbs: starting barbs {[x.name for x in self.data['barb']]}")
        # skip interval intended for 2-D fields
        if len(self.x.shape) == 2:
            cs2 = self.m.barbs(self.x[::skip,::skip], self.y[::skip,::skip], self.data['barb'][0][::skip,::skip], self.data['barb'][1][::skip,::skip],
                    alpha=alpha, length=5.5, linewidth=0.25, sizes={'emptybarb':0.05}, ax=self.ax)
        if len(self.x.shape) == 1:
            u2d = self.latlonGrid(self.data['barb'][0])
            v2d = self.latlonGrid(self.data['barb'][1])
            # rotate vectors so they represent the direction properly on the map projection
            u10_rot, v10_rot, x, y = self.m.rotate_vector(u2d, v2d, self.lon2d, self.lat2d, returnxy=True)
            cs2 = self.m.barbs(x[::skip,::skip], y[::skip,::skip], u10_rot[::skip,::skip], v10_rot[::skip,::skip], 
                    alpha=alpha, length=5.5, linewidth=0.25, sizes={'emptybarb':0.05}, ax=self.ax)
    
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
       data = self.data['contour'][0]
       data.values = ndimage.gaussian_filter(data, sigma=[0,4,4])
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

               logging.debug(f'member {member} creating axes at {x},{y}')
               thisax = fig.add_axes([x,y,w,h])

               thisax.axis('on')
               for axis in ['top','bottom','left','right']: thisax.spines[axis].set_linewidth(0.5)
               self.m.drawcoastlines(ax=thisax, linewidth=0.3)
               self.m.drawstates(linewidth=0.15, ax=thisax)
               self.m.drawcountries(ax=thisax, linewidth=0.3)
               thisax.text(0.03,0.97,member+1,ha="left",va="top",bbox=dict(boxstyle="square",lw=0.5,fc="white"), transform=thisax.transAxes)
               
               # plot, unless file that has fill field is missing, then skip
               if member not in self.missing_members and member < self.ENS_SIZE:
                   data = self.latlonGrid(self.data['fill'][0].isel(mem=memberidx))
                   logging.debug(f"plotStamp: starting contourf with regridded array member {memberidx}")
                   cs1 = self.m.contourf(self.x2d, self.y2d, data, levels=levels, cmap=cmap, norm=norm, extend='max', ax=thisax)
                   memberidx += 1

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
       fig.text(0.51, 0.22 - 25/float(fig_height_px), initstr, transform=fig.transFigure)
       fig.text(0.51, 0.22 - 40/float(fig_height_px), validstr, transform=fig.transFigure)

       # add NCAR logo and text below logo
       x, y = fig.transFigure.transform((0.51,0))
       fig.figimage(plt.imread('ncar.png'), xo=x, yo=y+15, zorder=1000)
       #plt.text(x+10, y+5, 'ensemble.ucar.edu', fontdict={'size':9, 'color':'#505050'}, transform=None)


    def getInitValidStr(self):
       initstr  = self.initdate.strftime(' Init: %a %Y-%m-%d %H UTC')
       shr = min(self.fhr)
       ehr = max(self.fhr)
       fmt = '%a %Y-%m-%d %H UTC'
       validstr = "Valid: " + (self.initdate+datetime.timedelta(hours=shr)).strftime(fmt)
       if ehr > shr: # range of valid times
            validstr += " - " + (self.initdate+datetime.timedelta(hours=ehr)).strftime(fmt)
       return initstr, validstr

    def saveFigure(self, trans=False):
        # place NCAR logo 57 pixels below bottom of map, then save image 
        if 'ensprod' in self.opts['fill']:  # CSS needed incase not a fill plot
           if not trans and not self.opts['fill']['ensprod'].endswith('stamp'):
             x, y = self.ax.transAxes.transform((0,0))
             self.fig.figimage(plt.imread('ncar.png'), xo=x, yo=(y-44))
             #plt.text(x+10, y-54, 'ensemble.ucar.edu', fontdict={'size':9, 'color':'#505050'}, transform=None)

        plt.savefig(self.outfile, dpi=90, transparent=trans)
        
        if self.opts['convert']:
            #command = 'convert -colors 255 %s %s'%(self.outfile, self.outfile)
            if not self.opts['fill']: ncolors = 48 #if no fill field exists
            elif self.opts['fill']['ensprod'] in ['prob', 'neprob', 'probgt', 'problt', 'neprobgt', 'neproblt']: ncolors = 48
            elif self.opts['fill']['name'] in ['crefuh']: ncolors = 48
            else: ncolors = 255
            #command = '/glade/u/home/sobash/pngquant/pngquant %d %s --ext=.png --force'%(ncolors,self.outfile)
            if os.environ['NCAR_HOST'] == "cheyenne":
                quant = '/glade/u/home/ahijevyc/bin_cheyenne/pngquant'
            else:
                quant = '/glade/u/home/ahijevyc/bin/pngquant'
            command = f"{quant} {ncolors} {self.outfile} --ext=.png --force"
            ret = subprocess.check_call(command.split())
        plt.clf()
        logging.info(f"created {self.outfile}")

def parseargs():
    '''Parse arguments and return dictionary of fill, contour and barb field parameters'''

    parser = argparse.ArgumentParser(description='Web plotting script for NCAR ensemble')
    parser.add_argument('date', help='initialization datetime')
    parser.add_argument('--autolevels', action='store_true', help='use min/max to determine levels for plot')
    parser.add_argument('-b', '--barb', help='barb field (FIELD_PRODUCT_THRESH)')
    parser.add_argument('-c', '--contour', help='contour field (FIELD_PRODUCT_THRESH)')
    parser.add_argument('-con', '--convert', default=True, action='store_false', help='run final image through imagemagick')
    parser.add_argument('-d', '--debug', action='store_true', help='turn on debugging')
    parser.add_argument('-f', '--fill', help='fill field (FIELD_PRODUCT_THRESH), field keys:'+','.join(list(fieldinfo.keys())))
    parser.add_argument('--fhr', nargs='+', type=float, default=[12], help='list of forecast hours')
    parser.add_argument('--meshstr', type=str, default='uni', help='mesh id or path to defining mesh')
    parser.add_argument('--nbarbs', type=int, default=50, help='max barbs in one dimension')
    parser.add_argument('--nlon_max', type=int, default=1500, help='max pts in longitude dimension')
    parser.add_argument('--nlat_max', type=int, default=1500, help='max pts in latitude dimension')
    parser.add_argument('--over', default=False, action='store_true', help='plot as overlay (no lines, transparent, no convert)')
    parser.add_argument('-sig', '--sigma', default=2, help='smooth probabilities using gaussian smoother')
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
        fieldtype = fields[f]['ensprod']
        fieldname = fields[f]['name']
        if fieldtype in ['prob', 'neprob', 'problt', 'probgt', 'neprobgt', 'neproblt', 'prob3d']: thresh = fields[f]['thresh']
        if fieldtype.startswith('mem'): member = int(fieldtype[3:])
        
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

            logging.info(f"perform {fieldtype} on {data.shape} data")
            if (fieldtype == 'mean'):
                data = data.mean(dim=["Time","mem"], keep_attrs=True)
            elif (fieldtype == 'pmm'): 
                data = compute_pmm(data)
            elif (fieldtype == 'max'):
                data = data.max(dim=["Time","mem"], keep_attrs=True)
            elif (fieldtype == 'min'):
                data = data.min(dim=["Time","mem"], keep_attrs=True)
            elif (fieldtype == 'var'):
                data = data.std(dim=["Time","mem"], keep_attrs=True)
            elif (fieldtype == 'maxstamp'):
                data = data.max(dim="Time", keep_attrs=True)
            elif (fieldtype == 'meanstamp'):
                data = data.mean(dim="Time", keep_attrs=True)
            elif (fieldtype == 'summean'):
                data = data.sum(dim="Time", keep_attrs=True)
                data = data.mean(dim="mem", keep_attrs=True)
            elif (fieldtype == 'maxmean'):
                data = data.max(dim="Time", keep_attrs=True)
                data = data.mean(dim="mem", keep_attrs=True)
            elif (fieldtype == 'summax'):
                data = data.sum(dim="Time", keep_attrs=True)
                data = data.max(dim="mem", keep_attrs=True)
            elif (fieldtype[0:3] == 'mem'):
                data = data.sel(mem=member)
            elif 'prob' in fieldtype:
                if fieldtype.endswith('prob') or fieldtype.endswith('gt'): data = (data>=thresh).astype('float')
                elif fieldtype.endswith('lt'): data = (data<thresh).astype('float')
            elif (fieldtype in ['prob3d']):
                data = (data>=thresh).astype('float')
                data = compute_prob3d(data, roi=14, sigma=float(fields['sigma']), type='gaussian')
            logging.debug(f'field {fieldname} has shape {data.shape} range {data.min()}-{data.min()}')

            datadict[f].append(data)

    return datadict, missing_list





def saveNewMap(plot, pk_file):
    logging.info(f"saveNewMap: pk_file={pk_file}")
    ll_lat, ll_lon, ur_lat, ur_lon = domains[plot.domain]['corners']
    lat_1, lat_2, lon_0 = 32., 46., -101.
    fig_width = 1080
    dpi = 90 
    fig = plt.figure(dpi=dpi)
    m = Basemap(projection='lcc', resolution='i', llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat, \
                lat_1=lat_1, lat_2=lat_2, lon_0=lon_0, area_thresh=1000)

    # compute height based on figure width, map aspect ratio, then add some vertical space for labels/colorbar
    fig_width  = fig_width/float(dpi)
    fig_height = fig_width*m.aspect + 0.93
    figsize = (fig_width, fig_height)
    fig.set_size_inches(figsize)

    # place map 0.7" from bottom of figure, leave rest of 0.93" at top for title (needs to be in figure-relative coords)
    x,y,w,h = 0.01, 0.7/float(fig_height), 0.98, 0.98*fig_width*m.aspect/float(fig_height)
    ax = fig.add_axes([x,y,w,h])
    for i in list(ax.spines.values()): i.set_linewidth(0.5)

    m.drawcoastlines(linewidth=0.5, ax=ax)
    m.drawstates(linewidth=0.25, ax=ax)
    m.drawcountries(ax=ax)
    m.drawcounties(linewidth=0.1, color='gray', ax=ax)

    # lat/lons from mpas_mesh file
    mpas_mesh = plot.mpas_mesh
    # min_grid_spacing_km used for grid spacing of interpolated lat-lon grid.
    min_grid_spacing_km = 2. * np.sqrt(plot.mpas_mesh["areaCell"].min()/np.pi)/1000

    lats = mpas_mesh["latCell"]
    lons = mpas_mesh["lonCell"]
    delta_deg = min_grid_spacing_km / 111
    if lons.max() > 180:
        lons[lons>180] -= 360
    # Replace m.lonmax with ur_lon? Otherwise, risk getting +179.999 when NW corner crosses the datetime
    nlon = int((m.lonmax - m.lonmin)/delta_deg)
    nlat = int((m.latmax - m.latmin)/delta_deg)
    nlon = np.clip(nlon, 1, plot.nlon_max)
    nlat = np.clip(nlat, 1, plot.nlat_max)
    lon2d, lat2d = np.meshgrid(np.linspace(m.lonmin, m.lonmax, nlon), np.linspace(m.latmin,m.latmax,nlat))
    # Convert to map coordinates instead of latlon to avoid latlon=True in contour and barb methods.
    x2d, y2d = m(lon2d,lat2d)
    # ibox: subscripts within lat/lon box. speed up triangulation in interp_weights
    ibox = (m.lonmin-1 <= lons ) & (lons < m.lonmax+1) & (m.latmin-1 <= lats) & (lats < m.latmax+1)
    lons = lons[ibox]
    lats = lats[ibox]
    x, y = m(lons,lats)

    logging.debug(f"saveNewMap: triangulate {len(lons)} pts")
    vtx, wts = interp_weights(np.vstack((lons,lats)).T,np.vstack((lon2d.flatten(), lat2d.flatten())).T)

    pickle.dump((fig,ax,m,lons,lats,delta_deg,lon2d,lat2d,x2d,y2d,ibox,x,y,vtx,wts), open(pk_file, 'wb'))
    


def drawOverlay(plot, pk_file):
    logging.info(f"drawOverlay: pk_file={pk_file}")
    ll_lat, ll_lon, ur_lat, ur_lon = domains[plot.domain]['corners']
    lat_1, lat_2, lon_0 = 32., 46., -101.
    fig_width = 1080
    dpi = 90 
    fig = plt.figure(dpi=dpi)
    m = Basemap(projection='lcc', resolution='i', llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat, \
                lat_1=lat_1, lat_2=lat_2, lon_0=lon_0, area_thresh=1000)

    # compute height based on figure width, map aspect ratio, then add some vertical space for labels/colorbar
    fig_width  = fig_width/float(dpi)
    fig_height = fig_width*m.aspect + 0.93
    figsize = (fig_width, fig_height)
    fig.set_size_inches(figsize)

    # place map 0.7" from bottom of figure, leave rest of 0.93" at top for title (needs to be in figure-relative coords)
    x,y,w,h = 0.01, 0.7/float(fig_height), 0.98, 0.98*fig_width*m.aspect/float(fig_height)
    ax = fig.add_axes([x,y,w,h])

    #drawcounties doesnt work when called by itself? so have to drawcoastines first with lw=0
    m.drawcoastlines(linewidth=0, ax=ax)
    m.drawcounties(ax=ax)
    ax.axis('off')
    plt.savefig('overlay_counties_%s.png'%domain, dpi=90, transparent=True)

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
    if len(ensemble.shape) < 3:
        logging.error('compute_neprob: needs ensemble of 2D arrays, not 1D arrays')
        sys.exit(1)
    y,x = np.ogrid[-roi:roi+1, -roi:roi+1]
    kernel = x**2 + y**2 <= roi**2
    ens_roi = ndimage.filters.maximum_filter(ensemble, footprint=kernel[np.newaxis,:])

    ens_mean = np.nanmean(ens_roi, axis=0)
    #ens_mean = np.nanmean(ensemble, axis=0)

    if type == 'uniform':
        y,x = np.ogrid[-sigma:sigma+1, -sigma:sigma+1]
        kernel = x**2 + y**2 <= sigma**2
        ens_mean = ndimage.filters.convolve(ens_mean, kernel/float(kernel.sum()))
    elif type == 'gaussian':
        ens_mean = ndimage.filters.gaussian_filter(ens_mean, sigma)
    else:
        logging.error(f"compute_neprob: unknown filter {type}")
        sys.exit(1)
    return ens_mean

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
