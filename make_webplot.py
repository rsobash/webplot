import logging
import os
import pdb
import sys
import time
from webplot import webPlot

logging.info('Begin Script')
stime = time.time()

regions = ['CONUS', 'NGP', 'SGP', 'CGP', 'MATL', 'NE', 'NW', 'SE', 'SW']



for dom in regions:
    newPlot = webPlot(domain=dom)
    logging.info('Reading Data')
    newPlot.readEnsemble()
    newPlot.loadMap()

    logging.info('Plotting Data')
    if newPlot.opts['interp']:
      newPlot.plotInterp()
    else:
      newPlot.plotFields()
      newPlot.plotTitleTimes()

    logging.info('Writing Image')
    newPlot.saveFigure(trans=newPlot.opts['over'])

etime = time.time()
logging.info(f'End Plotting (took {etime-stime:.2f} sec)')

