import argparse
import logging
import numpy as np
import pandas as pd
import pdb
import pickle
import pygrib
import os
from scipy.spatial import KDTree, Delaunay
import xarray


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

class regrid:
    '''regrid MPAS'''
    def __init__(self):
        self.args = parseargs()

        # Source mesh
        path = self.args.init_file
        logging.info(f"get mpas mesh file {path}")
        mpas_mesh = xarray.open_dataset(path)
        lonCell = mpas_mesh['lonCell']
        lonCell = np.degrees(lonCell) #convert radians to degrees
        lonCell[lonCell>=180] -= 360
        mpas_mesh["lonCell"] = lonCell
        mpas_mesh["latCell"] = np.degrees(mpas_mesh["latCell"]) #convert radians to degrees

        # Source data
        fhr = self.args.fhr[0]
        idate = pd.to_datetime(self.args.date)
        valid_time = idate + pd.Timedelta(fhr, unit="hour")
        ifile = os.path.join("/glade/campaign/mmm/parc/schwartz/MPAS_JEDI/15-3km_mesh/cold_start",
                idate.strftime("%Y%m%d%H"),
                f'diag.{valid_time.strftime("%Y-%m-%d_%H.%M.%S")}.nc' )
        var = "refl10cm_1km_max"
        logging.info(f"open {ifile} {var}")
        data = xarray.open_dataset(ifile)[var]

        # grib file with destination mesh
        today = pd.Timestamp.now().strftime('%Y%m%d00') 
        grb_with_dest_mesh = f'/glade/scratch/sobash/HRRR/{today}/hrrr.t00z.wrfsfcf00.grib2'
        with pygrib.open(grb_with_dest_mesh) as f:
            grb = f.readline()
            dest_latlons = grb.latlons()
        # destination lats and lons
        dest_lats, dest_lons = dest_latlons
        # weights and vertices for regridding 
        pk_file = os.path.join(os.getenv('TMPDIR'), f"{self.args.meshstr}.pk")
        if not os.path.exists(pk_file):
            saveNewMap(mpas_mesh, dest_latlons, pk_file)
        logging.info(f"load {pk_file}")
        (ibox, vtx, wts) = pickle.load(open(pk_file, 'rb'))

        # drop cells far away from destination mesh (ibox=False)
        data = data.isel(Time=0,nCells=ibox)
        logging.info("interpolatetri(vtx and wts)")
        out = interpolatetri(data.values, vtx, wts)
        logging.info("reshape")
        out = np.reshape(out, dest_lats.shape)
        logging.info("DataArray")
        out = xarray.DataArray(data=out, name=var, 
                coords = dict(lat=(["y","x"], dest_lats), lon=(["y","x"], dest_lons)), 
                dims=["y","x"], attrs=data.attrs)
        logging.info("to_netcdf")
        out.to_netcdf('out.nc')

        ogrb = 'out.grb2'
        logging.info(f"write {ogrb}")
        grb.values = out.values
        grb.dataDate = int(idate.strftime("%Y%m%d"))
        grb.dataTime = int(idate.strftime("%H%M"))
        grb["forecastTime"] = fhr
        with open(ogrb, 'wb') as f:
            f.write(grb.tostring())

def interpolatetri(values, vtx, wts):
    return np.einsum('nj,nj->n', np.take(values, vtx), wts)

def parseargs():
    '''Parse arguments and return dictionary of parameters'''

    parser = argparse.ArgumentParser(description='regrid MPAS',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('date', help='model initialization time')
    parser.add_argument('init_file', help='path to file with latCell, lonCell of mesh')
    parser.add_argument('-d', '--debug', action='store_true', help='turn on debugging')
    parser.add_argument('--ENS_SIZE', type=int, default=1, help='number of members in ensemble')
    parser.add_argument('--fhr', nargs='+', type=float, default=[12], help='list of forecast hours')
    parser.add_argument('--idir', type=str, default='/glade/campaign/mmm/parc/schwartz/MPAS_ensemble_paper',
            help="path to model output")
    parser.add_argument('--meshstr', type=str, default="HRRR", help=f'mesh id. used to prefix pickle file')

    args = parser.parse_args()
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    assert os.path.exists(args.init_file), "--init_file must be path to file with latCell, lonCell, etc."

    return args

def saveNewMap(mpas_mesh, dest_latlons, pk_file):
    logging.info(f"saveNewMap: pk_file={pk_file}")

    # destination lats and lons
    dest_lats, dest_lons = dest_latlons

    # Source lats and lons
    latCell = mpas_mesh.latCell
    lonCell = mpas_mesh.lonCell
    # To save time, triangulate near destination mesh
    # ibox: subscripts. speed up triangulation in interp_weights
    tree = KDTree(np.c_[dest_lons.ravel(), dest_lats.ravel()])
    dd, ii = tree.query(np.c_[lonCell, latCell], distance_upper_bound=0.25)
    ibox = np.isfinite(dd)
    latCell = latCell[ibox]
    lonCell = lonCell[ibox]

    logging.info(f"saveNewMap: triangulate {len(lonCell)} pts")
    vtx, wts = interp_weights(np.vstack((lonCell,latCell)).T,np.vstack((dest_lons.flatten(), dest_lats.flatten())).T)

    pickle.dump((ibox,vtx,wts), open(pk_file, 'wb'))
    
# from https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids
def interp_weights(xyz, uvw):
    tri = Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    d = xyz.shape[1]
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

if __name__ == "__main__":
    regrid()
