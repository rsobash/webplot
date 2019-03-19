import matplotlib.pyplot as plt
from subprocess import call # to use "mogrify"
import datetime as dt
def mysavfig(ofile, string="", timestamp=True, size=5, **kwargs):
    if timestamp:
        string = string + '\n' + 'created '+str(dt.datetime.now(tz=None)).split('.')[0] # + '\n' # extra newline to keep timestamp onscreen.
    th = plt.annotate(string, xy=(2,1), xycoords='figure pixels', horizontalalignment='left', verticalalignment='bottom', size=size)
    #plt.tight_layout() # Tried this but it screwed up SHARPY_skewt. Also tried bbox_inches='tight' below but not in that version of matplotlib.
    # Tried this but it made a large blank space on the left side of SHARPY_skewt.
    #    plt.savefig(ofile,dpi=dpi,bbox_inches='tight')
    plt.savefig(ofile, **kwargs)
    th.remove() # permits the figure to reused without overlaying multiple "created . . " dates.
    cmd = "mogrify +matte -type Palette -colors 255 " + ofile # prevents flickering when displaying on yellowstone.
    print("created", ofile)
    return call(cmd.split()) 

