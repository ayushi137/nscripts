import os
from astropy.io import fits
import numpy as np

directory = '../moddragonflydata/'
subdir = '/photometered/'
caldir = '/calframes/'

objectdates = {}
objectdates['PGM_1_2'] = ['2013-09-23','2013-09-24','2013-09-25','2013-09-26',
                          '2013-09-27']
objectdates['spi1_1'] = ['2014-10-31','2014-11-01','2014-11-19',
                         '2014-11-22','2014-11-23','2014-11-24']
keys = objectdates.keys()

mvs = {}
dts = {}
fts = {}

for key in keys:
    print 'OBJECT '+key
    dates = objectdates[key]
    for date in dates:
        print 'DATE '+date
        di = directory+key+subdir+date+'/'
        cdi = directory+key+caldir+date+'/'
        dfiles = os.listdir(di)
        cfiles = os.listdir(cdi)
        dfiles.sort()
        cfiles.sort()
        medianvals = {}
        datatimes = {}
        flattimes = {}
        for f in dfiles:
            if '_photo.fits' in f:
                data,dheader = fits.getdata(di+f,header = True)
                serialno = f[0:9]
                expn = f[10:12]
                ident = serialno+'_'+expn
                try:
                    datatimes[serialno].append(dheader['DATE'])
                except KeyError:
                    datatimes[serialno] = [dheader['DATE']]
                medianvals[ident] = np.median(data)
        for f in cfiles:
            if 'dark' not in f and 'DS' not in f:
                cheader = fits.getheader(cdi+f)
                serialno = f[11:20]
                try:
                    flattimes[serialno].append(cheader['DATE'])
                except KeyError:
                    flattimes[serialno] = [cheader['DATE']]
        mvs[date] = medianvals
        dts[date] = datatimes
        fts[date] = flattimes

