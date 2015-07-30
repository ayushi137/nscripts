import os
from astropy.io import fits
import numpy as np

directory = '../moddragonflydata/'
subdir = '/raw_lights/'

objectdates = {}
objectdates['PGM_1_2'] = ['2013-09-23','2013-09-24','2013-09-25','2013-09-26',
                          '2013-09-27']
objectdates['spi1_1'] = ['2014-10-31','2014-11-01','2014-11-19',
                         '2014-11-22','2014-11-23','2014-11-24']
keys = objectdates.keys()

for key in keys:
    print 'OBJECT '+key
    dates = objectdates[key]
    for date in dates:
        di = directory+key+subdir+date+'/'
        dfiles = os.listdir(di)
        dfiles.sort()
        for f in dfiles:
            if '.fits' in f:
                data,dheader = fits.getdata(di+f,header = True)
                serialno = f[0:9]
                try:
                    print date,serialno,dheader['FILTNAM']
                except KeyError:
                    print 'No filt ', date, serialno
