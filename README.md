# webplot plotting package

create graphics from MPAS ensemble

For example, to plot max precipitation accumulation with contours of 500hPa mean height, and 
mean 500hPa wind barbs from forecast hours [12, 18]:

```
python webplot.py 20170518 --fill precipacc/max --fhr 12 18 --contour hgt500/mean \
    --barb wind500/mean --title 'Max precip acc, mean 500hPa hgt [m] and wind barbs [kt]'
```

MPAS initialization time provided as first argument.

```
usage: webplot.py [-h] [--autolevels] [-b BARB] [-c CONTOUR] [-con] [-d]
                  [--domain {CONUS,NA,SGP,NGP,CGP,SW,NW,SE,NE,MATL}] [-f FILL] [--fhr FHR [FHR ...]]
                  [--idir IDIR] [--meshstr MESHSTR] [--nbarbs NBARBS] [--nlon_max NLON_MAX]
                  [--nlat_max NLAT_MAX] [--sigma SIGMA] [-t TITLE]
                  date

Web plotting script for NCAR ensemble

positional arguments:
  date                  model initialization time

options:
  -h, --help            show this help message and exit
  --autolevels          use min/max to determine levels for plot (default: False)
  -b BARB, --barb BARB  barb field (FIELD_PRODUCT_THRESH) (default: None)
  -c CONTOUR, --contour CONTOUR
                        contour field (FIELD_PRODUCT_THRESH) (default: None)
  -con, --convert       run final image through imagemagick (default: True)
  -d, --debug           turn on debugging (default: False)
  --domain {CONUS,NA,SGP,NGP,CGP,SW,NW,SE,NE,MATL}
                        domain of plot (default: CONUS)
  -f FILL, --fill FILL  fill field (FIELD_PRODUCT_THRESH), FIELD options:precip,precip-24hr,precip-
                        48hr,precipacc,sbcape,mlcape,mucape,sbcinh,mlcinh,pwat,t2,t2depart,t2-
                        0c,mslp,td2,td2depart,thetae,thetapv,rh2m,pblh,hmuh,hmneguh,hmuh03,hmuh01,rvort1,sspf,
                        hmup,hmdn,hmwind,hmgrp,cref,ref1km,srh3,srh1,shr06mag,shr01mag,zlfc,zlcl,ltg1,ltg2,ltg
                        3,olrtoa,thck1000-500,thck1000-
                        850,hgt200,hgt250,hgt300,hgt500,hgt700,hgt850,hgt925,speed200,speed250,speed300,speed5
                        00,speed700,speed850,speed925,temp200,temp250,temp300,temp500,temp700,temp850,temp925,
                        td500,td700,td850,td925,rh300,rh500,rh700,rh850,rh925,pvort320k,speed10m,speed10m-
                        tc,stp,uhratio,crefuh,wind10m,windsfc,wind1km,wind6km,windpv,shr06,shr01,wind200,wind2
                        50,wind300,wind500,wind700,wind850,wind925,vort500,vort700,vort850,vortpv PRODUCT may
                        be one of [max,maxstamp,min,mean,meanstamp,prob,neprob,problt,neproblt,paintball,stamp
                        ,spaghetti] (default: None)
  --fhr FHR [FHR ...]   list of forecast hours (default: [12])
  --idir IDIR           path to model output (default: /glade/campaign/mmm/parc/schwartz/MPAS_ensemble_paper)
  --meshstr MESHSTR     mesh id ['uni', '15-3km_mesh', '15km_mesh'] or path to file with lat/lon/area of mesh
                        (default: 15km_mesh)
  --nbarbs NBARBS       max barbs in one dimension (default: 32)
  --nlon_max NLON_MAX   max pts in longitude dimension (default: 1500)
  --nlat_max NLAT_MAX   max pts in latitude dimension (default: 1500)
  --sigma SIGMA         smooth probabilities using gaussian smoother (default: 2)
  -t TITLE, --title TITLE
                        title for plot (default: None)
```

### Installation

Use environment.yaml to create conda environment.

Use f2py3 to install vert2cell.

