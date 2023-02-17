import logging
import pdb
import time
from webplot import webPlot



stime = time.time()

regions = ['CONUS', 'NGP', 'SGP', 'CGP', 'MATL', 'NE', 'NW', 'SE', 'SW']
regions = ['CONUS', 'CGP']


for dom in regions:
    Plot = webPlot(domain=dom)
    logging.debug('Writing Image')
    Plot.saveFigure()

etime = time.time()
logging.info(f'End Plotting (took {etime-stime:.2f} sec)')

