import os
from numpy import *
from callastrometry import callastrometry, scrubwcsheader

directory = '../moddragonflydata/'
subdir = '/NCalibrated/'

objectdates = {}
objectdates['PGM_1_2'] = ['2013-09-23','2013-09-24','2013-09-25','2013-09-26',
                          '2013-09-27']
objectdates['spi1_1'] = ['2014-10-31','2014-11-01','2014-11-19',
                         '2014-11-22','2014-11-23','2014-11-24']

keys = objectdates.keys()

for key in keys:
    dates = objectdates[key]
    for date in dates:
        di = directory+key+subdir+date+'/'
        files = os.listdir(di)
        files.sort()
        for f in files:
            if '.fits' in f:
                print 'FILE ',di+f
                callastrometry(di+f)
