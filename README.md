# webplot plotting package

This is a simplified version of the webplot plotting library used to create graphics for the NCAR ensemble system.

For example, to create a 2-m ensemble mean temperature plot using ensemble data initialized at 2017061900, valid at 2017061912:

make_webplot.py -d=2017061900 -f=t2_mean -tr=12 -t='2-m Ensemble Mean Temperature (F)'
