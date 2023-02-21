# webplot plotting package

create graphics from MPAS ensemble

For example, to plot max precipitation accumulation, mean 500hPa vorticity, and mean wind barbs 
from the valid period [2017061912, 2017061918]:

```
python webplot.py 20170619 --fill precipacc/max --fhr 12 18 --contour vort500/mean --barb wind500/mean
    --title 'max precip accumuation, mean 500hPa vort, mean 500hPa wind barbs'
```

MPAS initialization time provided as first argument.

```
usage: webplot.py [-h] [--autolevels] [-b BARB] [-c CONTOUR] [-con] [-d]
                  [--domain {CONUS,NA,SGP,NGP,CGP,SW,NW,SE,NE,MATL}] [-f FILL] [--fhr FHR [FHR ...]]
                  [--meshstr MESHSTR] [--nbarbs NBARBS] [--nlon_max NLON_MAX] [--nlat_max NLAT_MAX]
                  [--sigma SIGMA] [-t TITLE]
                  date

Web plotting script for NCAR ensemble

positional arguments:
  date                  initialization datetime

options:
  -h, --help            show this help message and exit
  --autolevels          use min/max to determine levels for plot
  -b BARB, --barb BARB  barb field (FIELD_PRODUCT_THRESH)
  -c CONTOUR, --contour CONTOUR
                        contour field (FIELD_PRODUCT_THRESH)
  -con, --convert       run final image through imagemagick
  -d, --debug           turn on debugging
  --domain {CONUS,NA,SGP,NGP,CGP,SW,NW,SE,NE,MATL}
                        domain of plot
  -f FILL, --fill FILL  fill field (FIELD_PRODUCT_THRESH), FIELD options:precip,precip-24hr,precip-
                        48hr,precipacc,sbcape,mlcape,mucape,sbcinh,mlcinh,pwat,t2,t2depart,t2-
                        0c,mslp,td2,td2depart,thetae,thetapv,rh2m,pblh,hmuh,hmneguh,hmuh03,hmuh01,rvort1,sspf,hmu
                        p,hmdn,hmwind,hmgrp,cref,ref1km,srh3,srh1,shr06mag,shr01mag,zlfc,zlcl,ltg1,ltg2,ltg3,olrt
                        oa,thck1000-500,thck1000-
                        850,hgt200,hgt250,hgt300,hgt500,hgt700,hgt850,hgt925,speed200,speed250,speed300,speed500,
                        speed700,speed850,speed925,temp200,temp250,temp300,temp500,temp700,temp850,temp925,td500,
                        td700,td850,td925,rh300,rh500,rh700,rh850,rh925,pvort320k,speed10m,speed10m-
                        tc,stp,uhratio,crefuh,wind10m,windsfc,wind1km,wind6km,windpv,shr06,shr01,wind200,wind250,
                        wind300,wind500,wind700,wind850,wind925,vort500,vort700,vort850,vortpv PRODUCT may be one
                        of
                        [max,maxstamp,min,mean,meanstamp,prob,neprob,problt,neproblt,paintball,stamp,spaghetti]
  --fhr FHR [FHR ...]   list of forecast hours
  --meshstr MESHSTR     mesh id or path to defining mesh
  --nbarbs NBARBS       max barbs in one dimension
  --nlon_max NLON_MAX   max pts in longitude dimension
  --nlat_max NLAT_MAX   max pts in latitude dimension
  --sigma SIGMA         smooth probabilities using gaussian smoother
  -t TITLE, --title TITLE
                        title for plot
```
