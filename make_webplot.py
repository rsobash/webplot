import logging
import os
import pdb
import sys
import time
from webplot import webPlot



stime = time.time()

regions = ['CONUS', 'NGP', 'SGP', 'CGP', 'MATL', 'NE', 'NW', 'SE', 'SW']
regions = ['CONUS', 'NGP']


for dom in regions:
    Plot = webPlot(domain=dom)
    logging.debug('Writing Image')
    Plot.saveFigure(trans=Plot.opts['over'])

etime = time.time()
logging.info(f'End Plotting (took {etime-stime:.2f} sec)')

