import numpy as np

def confidence_interval(x):
    nsz = x.size
    std = np.std(x)
    # From 11 Aug 2016 email from Fanglin
    # 95% confidence level defined by {-intvl, intvl}
    if nsz>=80:
        intvl=1.960*std/np.sqrt(nsz-1)
    if nsz>=40 and nsz <80:
        intvl=2.000*std/np.sqrt(nsz-1)
    if nsz>=20 and nsz <40:
        intvl=2.042*std/np.sqrt(nsz-1)
    if nsz<20:
        intvl=2.228*std/np.sqrt(nsz-1)

    return intvl
    

