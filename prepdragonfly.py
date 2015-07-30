from astropy.io import fits
import os
from numpy import *
from callastrometry import scrubwcsheader

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
            print 'FILE ', di+f
            scrubwcsheader(di+f)
            i = fits.open(di+f)
            idata = i[0].data
            ihder = i[0].header
            i.close()
            ihder['UNITS'] = 'ADU'
            try:
                if ihder['FLIPLR'] == 'RA_RL':
                    continue
                elif ihder['FLIPLR'] == 'RA_LR':
                    idata = fliplr(idata)
                    ihder['FLIPLR'] = 'RA_RL'
                    fits.writeto(di+f,idata,ihder,clobber = True)
            except KeyError:
                idata = fliplr(idata)
                ihder['FLIPLR'] = 'RA_RL'
                fits.writeto(di+f,idata,ihder,clobber = True)
