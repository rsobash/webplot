from fieldinfo import fieldinfo,readNCLcm
import numpy as np

# fieldinfo should have been imported from fieldinfo module.
# Copy fieldinfo dictionary for MPAS. Change some fnames.
fieldinfo['td2']['fname'] = ['surface_dewpoint']
fieldinfo['td2depart']['fname'] = ['surface_dewpoint']
for plev in ['500', '700', '850']:
    fieldinfo['vort'+plev] = {'levels' : np.array([0,9,12,15,18,21,24,27,30,33])*1e-5, 'cmap': readNCLcm('prcp_1'), 'fname': ['vorticity_'+plev+'hPa']}
fieldinfo['vortpv']        = {'levels' : [0,9,12,15,18,21,24,27,30,33], 'cmap': readNCLcm('prcp_1'), 'fname': ['vort_pv']}
fieldinfo['stp']['fname'] = ['sbcape','lcl','srh_0_1km','uzonal_6km','umeridional_6km','uzonal_surface','umeridional_surface']
