import sys
import pdb
#sys.path.append('/glade/u/home/ahijevyc/lib/python2.7/site-packages/SHARPpy-1.3.0-py2.7.egg')
#sys.path.append('/glade/u/home/ahijevyc/lib/python2.7/site-packages')
import sharppy
import sharppy.sharptab.interp as interp
import sharppy.sharptab.params as params
import sharppy.sharptab.profile as profile
import sharppy.sharptab.thermo as thermo
import sharppy.sharptab.utils as utils
import sharppy.sharptab.winds as winds

import numpy as np
from io import StringIO
import matplotlib.pyplot as plt
import cartopy


def parseCLS(sfile):
    ## read in the file
    data = np.array([l.strip() for l in sfile.split('\n')])

    latitude = float(data[3].split(',')[5])
    longitude = float(data[3].split(',')[4])
    ## necessary index points
    start_idx = 15
    finish_idx = len(data)
    
    ## put it all together for StringIO
    full_data = '\n'.join(data[start_idx : finish_idx][:])
    sound_data = StringIO( full_data )

    ## read the data into arrays
    data = np.genfromtxt( sound_data )
    clean_data = []
    for i in data:
        if i[1] != 9999 and i[2] != 999 and i[3] != 999 and i[7] != 999 and i[8] != 999 and i[14] != 99999:
            clean_data.append(i)
    p = np.array([i[1] for i in clean_data])
    h = np.array([i[14] for i in clean_data])
    T = np.array([i[2] for i in clean_data])
    Td = np.array([i[3] for i in clean_data])
    wdir = np.array([i[8] for i in clean_data])
    wspd = np.array([i[7] for i in clean_data])
    wspd = utils.MS2KTS(wspd)
    wdir[wdir == 360] = 0. # Some wind directions are 360. Like in /glade/p/work/ahijevyc/GFS/Joaquin/g132325165.frd

    max_points = 250
    s = p.size/max_points
    if s == 0: s = 1
    print("stride=",s)
    return p[::s], h[::s], T[::s], Td[::s], wdir[::s], wspd[::s], latitude, longitude

def thetas(tempC, presvals):
    return ((tempC + thermo.ZEROCNK) / (np.power((1000. / presvals),thermo.ROCP))) - thermo.ZEROCNK

	 
	  
def dewpoint_approximation(T,RH):
    # approximation valid for
    # 0 degC < T < 60 degC
    # 1% < RH < 100%
    # 0 degC < Td < 50 degC 
    # constants
    a = 17.271
    b = 237.7 # degC
	 
    Td = (b * gamma(T,RH)) / (a - gamma(T,RH))
 
    return Td
 
 
def gamma(T,RH):
    # constants
    a = 17.271
    b = 237.7 # degC
    
    g = (a * T / (b + T)) + np.log(RH/100.0)
 
    return g

def draw_background(ax, dry=np.arange(-50,110,20), moist=np.arange(-5,40,5), presvals=np.arange(1050, 0, -10), mix=[1,2,4,8,12,16,20,24]):

    # Plot dry adiabats
    for t in dry:
        ax.semilogy(thetas(t, presvals), presvals,'r-', alpha=0.2, linewidth=1)
    
    # Plot moist adiabats to topp
    topp = 200
    moist_presvals = np.arange(np.max(presvals), topp, -10)
    for t in moist:
        tw = []
        for p in moist_presvals:
            tw.append(thermo.wetlift(1000., t, p))
        ax.semilogy(tw, moist_presvals, 'k-', alpha=0.2, linewidth=0.5, linestyle='dashed')
        # add moist adiabat text labels
        ax.text(thermo.wetlift(1000., t, topp+15), topp+15, t, va='center', ha='center', size=5.6, alpha=0.3, clip_on=True)


    # Mixing ratio lines and labels to topp

    topp = 650
    for w in mix:
        ww = []
        for p in presvals:
            ww.append(thermo.temp_at_mixrat(w,p))
        ax.semilogy(ww, presvals, 'g', alpha=0.35, linewidth=1, linestyle="dotted")
        ax.text(thermo.temp_at_mixrat(w, topp), topp, w, va='bottom', ha='center', size=6.7, color='g', alpha=0.35)

    # Disables the log-formatting (10 to the power of x) that comes with semilogy
    ax.yaxis.set_major_formatter(plt.ScalarFormatter())
    ax.set_yticks(np.linspace(100,1000,10))
    ax.set_ylim(1050,100)
    plt.ylabel('Pressure (hPa)')

    # The first time this axis object is returned, no xtick gridlines show up for -110 to -60C
    ax.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax.set_xlim(-50,45)
    plt.xlabel('Temperature (C)')

    ax.grid(True, linestyle='solid', alpha=0.5)


def draw_hodo():
    bbox_props = dict(boxstyle="square", color="w", alpha=0.5, lw=0.5)
    ax = plt.axes([.5,.63,.2,.26])
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    az = np.radians(-35)
    for i in range(10,100,10):
        # Draw the range rings around the hodograph.
        lw = .4 if i % 20 == 0 else 0.25
        circle = plt.Circle((0,0),i,color='k', alpha=.3, fill=False, lw=lw)
        ax.add_artist(circle)
        if i % 20 == 0 and i < 100:
            # The minus 2 nudges it a little closer to origin. 
            plt.text((i-2)*np.cos(az),(i-2)*np.sin(az),str(i)+" kt",rotation=np.degrees(az),size=5,alpha=.4,
                    ha='center', zorder=1, bbox=bbox_props)

    ax.set_xlim(-40,80)
    ax.set_ylim(-60,60)
    ax.axhline(y=0, color='k')
    ax.axvline(x=0, color='k')

    return ax


def add_hodo(ax, prof, lw=1, color='black', ls='solid', size=5, AGLs=[1,2,6,10]):
    # Return 2D hodograph line and list of text AGL labels.


    # Draw hodo line
    u_prof = prof.u[prof.pres >= 100]
    v_prof = prof.v[prof.pres >= 100]
    hodo, = ax.plot(u_prof,v_prof, 'k-', lw=lw, color=color, ls=ls)

    # position AGL text labels
    bbox_props = dict(boxstyle="square", fc="w", ec="0.5", linewidth=0.5, alpha=0.7)
    AGL_labels = []
    AGL_labels.append(ax.text(-35,-55,'km AGL',size=size,bbox=bbox_props))
    for a in AGLs:
        in_meters = a*1000
        if in_meters <= np.max(interp.to_agl(prof, prof.hght)): 
            junk = ax.text(0,0,str(a),ha='center',va='center',size=size,bbox=bbox_props,color=color)
            ind = np.min(np.where(interp.to_agl(prof,prof.hght)>in_meters))
            junk.set_position((prof.u[ind],prof.v[ind]))
            junk.set_clip_on(True)
            AGL_labels.append(junk)
    return hodo, AGL_labels
   

def wind_barb_spaced(ax, prof, xpos=1.0, yspace=0.04):
    # yspace = 0.04 means ~26 barbs.
    s = []
    bot = 2000.
    # Space out wind barbs evenly on log axis.
    for ind, i in enumerate(prof.pres):
        if i < 100: break
        if np.log(bot/i) > yspace:
            s.append(ind)
            bot = i
    # x coordinate in (0-1); y coordinate in pressure log(p)
    b = plt.barbs(xpos*np.ones(len(prof.pres[s])), prof.pres[s], prof.u[s], prof.v[s],
          length=5.8, lw=0.35, clip_on=False, transform=ax.get_yaxis_transform())

    # 'knots' label under wind barb stack
    kts = ax.text(1.0, prof.pres[0]+10, 'knots', clip_on=False, transform=ax.get_yaxis_transform(),ha='center',va='top',size=7)
    return b, kts


def add_globe(longitude, latitude):
    # TODO: avoid matplotlib depreciation warning about creating a unique id for each axes instance. You need a new axes
    # instance whenever lat/lon changes.
    # Globe with dot on location.
    mapax = plt.axes([.795, 0.09,.18,.18], projection=cartopy.crs.Orthographic(longitude, latitude))
    mapax.add_feature(cartopy.feature.OCEAN, zorder=0)
    mapax.add_feature(cartopy.feature.LAND, zorder=0, linewidth=0) # linewidth=0 or coastlines are fuzzy
    mapax.set_global()
    sloc = mapax.plot(longitude, latitude,'o', color='green', markeredgewidth=0, markersize=4., transform=cartopy.crs.Geodetic())
    return mapax


def indices(prof, debug=False):

    # return a formatted-string list of stability and kinematic indices

    sfcpcl = params.parcelx(prof, flag=1)
    mupcl = params.parcelx(prof, flag=3) # most unstable
    mlpcl = params.parcelx(prof, flag=4) # 100 mb mean layer parcel


    pcl = mupcl
    sfc = prof.pres[prof.sfc]
    p3km = interp.pres(prof, interp.to_msl(prof, 3000.))
    p6km = interp.pres(prof, interp.to_msl(prof, 6000.))
    p1km = interp.pres(prof, interp.to_msl(prof, 1000.))
    mean_3km = winds.mean_wind(prof, pbot=sfc, ptop=p3km)
    sfc_6km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p6km)
    sfc_3km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p3km)
    sfc_1km_shear = winds.wind_shear(prof, pbot=sfc, ptop=p1km)
    #print "0-3 km Pressure-Weighted Mean Wind (kt):", utils.comp2vec(mean_3km[0], mean_3km[1])[1]
    #print "0-6 km Shear (kt):", utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1]
    srwind = params.bunkers_storm_motion(prof)
    srh3km = winds.helicity(prof, 0, 3000., stu = srwind[0], stv = srwind[1])
    srh1km = winds.helicity(prof, 0, 1000., stu = srwind[0], stv = srwind[1])
    #print "0-3 km Storm Relative Helicity [m2/s2]:",srh3km[0]

    #### Calculating variables based off of the effective inflow layer:

    # The effective inflow layer concept is used to obtain the layer of buoyant parcels that feed a storm's inflow.
    # Here are a few examples of how to compute variables that require the effective inflow layer in order to calculate them:

    stp_fixed = params.stp_fixed(sfcpcl.bplus, sfcpcl.lclhght, srh1km[0], utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1])
    ship = params.ship(prof)

    # If you get an error about not converting masked constant to python int 
    # use the round() function instead of int() - Ahijevych May 11 2016
    # 2nd element of list is the # of decimal places
    indices = {'SBCAPE': [sfcpcl.bplus, 0, 'J $\mathregular{kg^{-1}}$'],
           'SBCIN': [sfcpcl.bminus, 0, 'J $\mathregular{kg^{-1}}$'],
           'SBLCL': [sfcpcl.lclhght, 0, 'm AGL'],
           'SBLFC': [sfcpcl.lfchght, 0, 'm AGL'],
           'SBEL': [sfcpcl.elhght, 0, 'm AGL'],
           'SBLI': [sfcpcl.li5, 0, 'C'],
           'MLCAPE': [mlpcl.bplus, 0, 'J $\mathregular{kg^{-1}}$'],
           'MLCIN': [mlpcl.bminus, 0, 'J $\mathregular{kg^{-1}}$'],
           'MLLCL': [mlpcl.lclhght, 0, 'm AGL'],
           'MLLFC': [mlpcl.lfchght, 0, 'm AGL'],
           'MLEL': [mlpcl.elhght, 0, 'm AGL'],
           'MLLI': [mlpcl.li5, 0, 'C'],
           'MUCAPE': [mupcl.bplus, 0, 'J $\mathregular{kg^{-1}}$'],
           'MUCIN': [mupcl.bminus, 0, 'J $\mathregular{kg^{-1}}$'],
           'MULCL': [mupcl.lclhght, 0, 'm AGL'],
           'MULFC': [mupcl.lfchght, 0, 'm AGL'],
           'MUEL': [mupcl.elhght, 0, 'm AGL'],
           'MULI': [mupcl.li5, 0, 'C'],
           '0-1 km SRH': [srh1km[0], 0, '$\mathregular{m^{2}s^{-2}}$'],
           '0-1 km Shear': [utils.comp2vec(sfc_1km_shear[0], sfc_1km_shear[1])[1], 0, 'kt'],
           '0-3 km SRH': [srh3km[0], 0, '$\mathregular{m^{2}s^{-2}}$'],
           '0-6 km Shear': [utils.comp2vec(sfc_6km_shear[0], sfc_6km_shear[1])[1], 0, 'kt'],
           'PWV': [params.precip_water(prof), 2, 'inch'],
           'K-index': [params.k_index(prof), 0, ''],
           'STP(fix)': [stp_fixed, 1, ''],
           'SHIP': [ship, 1, '']}



    eff_inflow = params.effective_inflow_layer(prof)
    if any(eff_inflow):
        ebot_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
        etop_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[1]))
        #print "Effective Inflow Layer Bottom Height (m AGL):", ebot_hght
        #print "Effective Inflow Layer Top Height (m AGL):", etop_hght
        effective_srh = winds.helicity(prof, ebot_hght, etop_hght, stu = srwind[0], stv = srwind[1])
        indices['Eff. SRH'] = [effective_srh[0], 0, '$\mathregular{m^{2}s^{-2}}$']
        #print "Effective Inflow Layer SRH (m2/s2):", effective_srh[0]
        ebwd = winds.wind_shear(prof, pbot=eff_inflow[0], ptop=eff_inflow[1])
        ebwspd = utils.mag( *ebwd )
        indices['EBWD'] = [ebwspd,0, 'kt']
        #print "Effective Bulk Wind Difference:", ebwspd
        scp = params.scp(mupcl.bplus, effective_srh[0], ebwspd)
        indices['SCP'] = [scp, 1, '']
        stp_cin = params.stp_cin(mlpcl.bplus, effective_srh[0], ebwspd, mlpcl.lclhght, mlpcl.bminus)
        indices['STP(cin)'] = [stp_cin, 1, '']
        #print "Supercell Composite Parameter:", scp
        #print "Significant Tornado Parameter (w/CIN):", stp_cin
        #print "Significant Tornado Parameter (fixed):", stp_fixed

    # Update the indices within the indices dictionary on the side of the plot.
    string = ''
    for index, value in sorted(indices.items()):
        if np.ma.is_masked(value[0]):
            if debug:
                print("skipping masked value for index=",index)
            continue
        if debug:
            print("index=",index)
            print("value=",value)
        format = '%.'+str(value[1])+'f'
        string += index + ": " +  format % value[0] + " " + value[2] + '\n'

    return string

