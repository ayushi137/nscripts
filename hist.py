import os
from astropy.io import fits
import matplotlib.pyplot as plt
from numpy import *
plt.ion()

chosenkey = 'FWHM'
outputdir = {}
outputdir['M0']='m0stats/'
outputdir['FWHM']='fwhmstats/'
xlabels = {}
xlabels['M0'] = 'm0'
xlabels['FWHM'] = 'fwhm ["]'

directory = '/mnt/gfsproject/naiad/njones/moddragonflydata/'
subdir = '/NCalibrated/'
outdir = '/backgroundsub/'
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
            if '35.2' not in f and '17.6' not in f:
                try:
                    print 'FILE ',pdi+f
                    serialno = f[0:9]
                    expn = f[10:12]
                    header = fits.getheader(pdi+f)
                    colour[serialno] = header['FILTNAM']
                    try:
                        datedict[serialno+'_'+expn].append(float(header[chosenkey]))
                    except (KeyError,NameError):
                        datedict[serialno+'_'+expn] = [float(header[chosenkey])]
                    try:
                        camdict[serialno].append(float(header[chosenkey]))
                    except (KeyError,NameError):
                        camdict[serialno] = [float(header[chosenkey])]
                except KeyError:
                    continue
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
    red = [item for item in red if item < 20]
    green = [item for sublist in green for item in sublist]
    green = [item for item in green if item < 20]
    try:
        plt.figure()
        plt.hist(red,bins = nbins)
        plt.title(objects+' - red')
        plt.ylabel('Number')
        plt.xlabel(xlabels[chosenkey])
        plt.savefig(outputdir[chosenkey]+objects+'_red.png')
        plt.close()
        plt.figure()
        plt.hist(green,bins = nbins)
        plt.title(objects+' - green')
        plt.ylabel('Number')
        plt.xlabel(xlabels[chosenkey])
        plt.savefig(outputdir[chosenkey]+objects+'_green.png')
        plt.close()
    except ValueError:
        continue

RED = array([item for sublist in RED for item in sublist])
RED = array([item for item in RED if item < 20])
GREEN = array([item for sublist in GREEN for item in sublist])
GREEN = array([item for item in GREEN if item < 20])

plt.figure()
plt.hist(RED,bins = nbins)
plt.title('Red Cameras')
plt.ylabel('Number')
plt.xlabel(xlabels[chosenkey])
plt.savefig(outputdir[chosenkey]+'redcams.png')
plt.close()
plt.figure()
plt.hist(GREEN,bins = nbins)
plt.title('Green Cameras')
plt.ylabel('Number')
plt.xlabel(xlabels[chosenkey])
plt.savefig(outputdir[chosenkey]+'greencams.png')
plt.close()
