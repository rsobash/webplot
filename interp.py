import argparse
import logging
import numpy as np
import pandas as pd
import pdb
import pickle
import pygrib
import os
from scipy.spatial import KDTree, Delaunay
import sys
import xarray


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

class regrid:
    '''regrid MPAS'''
    def __init__(self):
        self.args = parseargs()
        idate = pd.to_datetime(self.args.date)
        idir = self.args.idir
        var = self.args.var

        # Source mesh
        path = self.args.init_file
        logging.info(f"get mpas mesh file {path}")
        mpas_mesh = xarray.open_dataset(path)
        lonCell = mpas_mesh['lonCell']
        lonCell = np.degrees(lonCell) #radians to degrees
        lonCell[lonCell>=180] -= 360
        mpas_mesh["lonCell"] = lonCell
        mpas_mesh["latCell"] = np.degrees(mpas_mesh["latCell"]) #radians to degrees

        grb_with_dest_mesh = self.args.grb_with_dest_mesh
        with pygrib.open(grb_with_dest_mesh) as f:
            grb = f.readline()
            dest_latlons = grb.latlons()
        grb.dataDate = int(idate.strftime("%Y%m%d"))
        grb.dataTime = int(idate.strftime("%H%M"))
        # destination lats and lons
        dest_lats, dest_lons = dest_latlons

        # weights and vertices for regridding 
        pk_file = os.path.join(os.getenv('TMPDIR'), f"{self.args.meshstr}.pk")
        if not os.path.exists(pk_file):
            saveNewMap(mpas_mesh, dest_latlons, pk_file)
        logging.info(f"load {pk_file}")
        (ibox, vtx, wts) = pickle.load(open(pk_file, 'rb'))

        # Source data
        valid_times = [idate + pd.Timedelta(fhr, unit="hour") for fhr in self.args.fhr]
        ifiles = [  os.path.join(idir,
                    idate.strftime("%Y%m%d%H"),
                    f'diag.{valid_time.strftime("%Y-%m-%d_%H.%M.%S")}.nc')
                    for valid_time in valid_times ]

        nt = len(ifiles) # used to reshape output
        logging.info(f"open {len(ifiles)} files")
        data = xarray.open_mfdataset(ifiles, concat_dim="Time", combine="nested")
        data = data.assign_coords(Time=valid_times)
        data = data[var].isel(nCells=ibox)
        # dot product of nCells x 3 arrays vtx and wtx, preserving time dimension
        out = np.einsum('tnj,nj->tn', np.take(data.values, vtx, axis=1), wts)
        out = np.reshape(out, (nt,*dest_lats.shape))
        out = xarray.DataArray(data=out, name=var, 
                coords = dict(
                    Time=data.Time,
                    lat=(["y","x"], dest_lats), 
                    lon=(["y","x"], dest_lons)), 
                dims=["Time","y","x"], attrs=data.attrs)

        ogrb = f'{idate.strftime("%Y%m%d_%H")}.grb2'
        logging.info(f"write {ogrb}")
        with open(ogrb, 'wb') as f:
            for i, fhr in enumerate(self.args.fhr):
                grb.values = out.values[i]
                grb["forecastTime"] = fhr
                f.write(grb.tostring())

        ocdf = f'{idate.strftime("%Y%m%d_%H")}.nc'
        logging.info(f"to_netcdf {ocdf}")
        xarray.concat(out, pd.Index(self.args.fhr, name="forecast_hour")).to_netcdf(ocdf)


def parseargs():
    '''Parse arguments and argparse Namespace'''

    parser = argparse.ArgumentParser(description='regrid MPAS',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('date', help='model initialization time')
    parser.add_argument('init_file', help='path to init.nc file with latCell, lonCell of mesh')
    parser.add_argument('-d', '--debug', action='store_true', help='turn on debugging')
    parser.add_argument('--ENS_SIZE', type=int, default=1, help='number of members in ensemble')
    parser.add_argument('--fhr', nargs='+', type=float, default=[12], help='list of forecast hours')
    today = pd.Timestamp.now().strftime('%Y%m%d00') 
    parser.add_argument('--grb_with_dest_mesh', default = f'/glade/scratch/sobash/HRRR/{today}/hrrr.t00z.wrfsfcf00.grib2',
            help = 'grib file with destination mesh')
    parser.add_argument('--idir', default='/glade/campaign/mmm/parc/schwartz/MPAS_ensemble_paper',
            help="path to model output")
    parser.add_argument('--meshstr', default="HRRR", help=f'mesh id. used to prefix pickle file')
    parser.add_argument('--var', default="refl10cm_1km_max", help=f'name of variable to regrid')

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
    assert wts.min() >= 0, 'got negative weight'

    pickle.dump((ibox,vtx,wts), open(pk_file, 'wb'))
    
# from https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids
def interp_weights(xyz, uvw):
    logging.info("interp_weights: Delaunay")
    tri = Delaunay(xyz)
    logging.info("interp_weights: find_simplex")
    simplex = tri.find_simplex(uvw)
    logging.info("interp_weights: vertices")
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    d = xyz.shape[1]
    delta = uvw - temp[:, d]
    logging.info("interp_weights: bary")
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

if __name__ == "__main__":
    regrid()
