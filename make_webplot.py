#!/usr/bin/env python

import sys, time, os
from webplot import webPlot, readGrid, drawOverlay, saveNewMap

def log(msg): print time.ctime(time.time()),':', msg

log('Begin Script'); stime = time.time()

regions = ['CONUS', 'NGP', 'SGP', 'CGP', 'MATL', 'NE', 'NW', 'SE', 'SW']

if not os.path.exists('picklefilename.pk'):
    saveNewMap(wrfout='wrfout_file_containing_lat_lons', domstr='name_for_domain') 

newPlot = webPlot()
log('Reading Data'); newPlot.readEnsemble()

for dom in regions:
    file_not_created, num_attempts = True, 0
    while file_not_created and num_attempts <= 3:
        newPlot.domain = dom

        newPlot.createFilename()
        fname = newPlot.outfile
 
        log('Loading Map for %s'%newPlot.domain)
        newPlot.loadMap()

        log('Plotting Data')
        if newPlot.opts['interp']:
          newPlot.plotInterp()
        else:
          newPlot.plotFields()
          newPlot.plotTitleTimes()

        log('Writing Image')
        newPlot.saveFigure(trans=newPlot.opts['over'])

        if os.path.exists(fname):
            file_not_created = False 
            log('Created %s, %.1f KB'%(fname,os.stat(fname).st_size/1000.0))
    
        num_attempts += 1

etime = time.time()
log('End Plotting (took %.2f sec)'%(etime-stime))

