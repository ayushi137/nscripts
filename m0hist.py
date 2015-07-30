import os
from astropy.io import fits
import matplotlib.pyplot as plt
from numpy import *
plt.ion()

outputdir = 'm0stats/'

directory = '/mnt/gfsproject/naiad/njones/moddragonflydata/'
subdir = '/NCalibrated/'
outdir = '/photometered/'
magdir = '/magnitudecalcs/'

objectdates = {}
objectdates['PGM_1_2'] = ['2013-09-23','2013-09-24','2013-09-25','2013-09-26',
                          '2013-09-27']
objectdates['spi1_1'] = ['2014-10-31','2014-11-01','2014-11-19',
                         '2014-11-22','2014-11-23','2014-11-24']

keys = objectdates.keys()
datedicts = {}
camdicts = {}
colours = {}
for key in keys:
    print 'OBJECT '+key
    dates = objectdates[key]
    camdict = {}
    colour = {}
    for date in dates:
        print 'DATE '+date
        di = directory+key+subdir+date+'/'
        pdi = directory+key+outdir+date+'/'
        mdi = directory+key+magdir+date+'/'
        Mdi = directory+key+magdir
        files = os.listdir(pdi)
        files.sort()
        datedict = {}
        for f in files:
            if 'convto' not in f and 'ADU' not in f:
                #print 'FILE ',pdi+f
                serialno = f[0:9]
                expn = f[10:12]
                header = fits.getheader(pdi+f)
                colour[serialno] = header['FILTNAM']
                try:
                    datedict[serialno+'_'+expn].append(float(header['M0']))
                except (KeyError,NameError):
                    datedict[serialno+'_'+expn] = [float(header['M0'])]
                try:
                    camdict[serialno].append(float(header['M0']))
                except (KeyError,NameError):
                    camdict[serialno] = [float(header['M0'])]
        datedicts[date] = datedict
    camdicts[key] = camdict
    colours[key] = colour
nbins = 10
RED = []
GREEN = []
for objects in camdicts.keys():
    red = []
    green = []
    for serialno in camdicts[objects].keys():
        if colours[objects][serialno] == 'SloanG':
            green.append(camdicts[objects][serialno])
            GREEN.append(camdicts[objects][serialno])
        elif colours[objects][serialno] == 'SloanR':
            red.append(camdicts[objects][serialno])
            RED.append(camdicts[objects][serialno])
        #try:
        #    plt.hist(camdicts[objects][serialno],bins = 5)
        #    plt.title(objects+'-'+serialno)
        #    plt.ylabel('Number')
        #    plt.xlabel('m_0')
        #    plt.figure()
        #except ValueError:
        #    'Single data point'
    red = [item for sublist in red for item in sublist]
    green = [item for sublist in green for item in sublist]
    plt.figure()
    plt.hist(red,bins = nbins)
    plt.title(objects+' - red')
    plt.ylabel('Number')
    plt.xlabel('m_0')
    plt.savefig(outputdir+objects+'_red.png')
    plt.figure()
    plt.hist(green,bins = nbins)
    plt.title(objects+' - green')
    plt.ylabel('Number')
    plt.xlabel('m_0')
    plt.savefig(outputdir+objects+'_green.png')

RED = array([item for sublist in RED for item in sublist])
GREEN = array([item for sublist in GREEN for item in sublist])

plt.figure()
plt.hist(RED,bins = nbins)
plt.title('Red Cameras')
plt.ylabel('Number')
plt.xlabel('m_0')
plt.savefig(outputdir+'redcams.png')
plt.figure()
plt.hist(GREEN,bins = nbins)
plt.title('Green Cameras')
plt.ylabel('Number')
plt.xlabel('m_0')
plt.savefig(outputdir+'greencams.png')
