import logging
import os
import pdb
import sys
import time
from webplot import webPlot, readGrid, saveNewMap

logging.info('Begin Script')
stime = time.time()

regions = ['CONUS', 'NGP', 'SGP', 'CGP', 'MATL', 'NE', 'NW', 'SE', 'SW']


newPlot = webPlot()
logging.info('Reading Data')
newPlot.readEnsemble()

for dom in regions:
    pk_file = newPlot.pk_file
    if not os.path.exists(pk_file):
        saveNewMap(newPlot, wrfout='latlonfile.nc', init_file='/glade/campaign/mmm/parc/ahijevyc/MPAS/uni/2018103000/init.nc')
    file_not_created, num_attempts = True, 0
    while file_not_created and num_attempts <= 3:
        newPlot.domain = dom

        newPlot.createFilename()
        fname = newPlot.outfile
 
        newPlot.loadMap()

        logging.info('Plotting Data')
        if newPlot.opts['interp']:
          newPlot.plotInterp()
        else:
          newPlot.plotFields()
          newPlot.plotTitleTimes()

        logging.info('Writing Image')
        newPlot.saveFigure(trans=newPlot.opts['over'])

        if os.path.exists(fname):
            file_not_created = False 
            logging.info(f'Created {fname} {os.stat(fname).st_size/1000:.1f} KB')
    
        num_attempts += 1

etime = time.time()
logging.info(f'End Plotting (took {etime-stime:.2f} sec)')

