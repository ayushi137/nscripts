import os
from astropy.io import fits
from astropy import wcs
import matplotlib.pyplot as plt
from numpy import *
from cartesian import cartesian
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord, AltAz
plt.ion()

def getAltAz(arr,header,time,location):
	soln = wcs.WCS(header)
	coords = cartesian([arange(arr.shape[1]),arange(arr.shape[0])])
	world = soln.wcs_pix2world(coords,0)
	radec = SkyCoord(ra=world[:,0],dec=world[:,1],frame='icrs',unit='deg')
	altaz = radec.transform_to(AltAz(obstime=time,location=telescope))
	return altaz.alt.deg,altaz.az.deg,coords[:,0],coords[:,1]

def airmass(alt):
    return (1./(sin(alt*(pi/180.))+0.50572*(6.07995+alt)**-1.6364))

lat = 32.902836
long = -105.528350
height = 2214
telescope = EarthLocation(lat=lat*u.deg,lon=long*u.deg,height = 2214*u.m)

chosenkey = 'FWHM'
outputdir = {}
outputdir['M0']='m0stats/'
outputdir['FWHM']='fwhmstats/'
xlabels = {}
xlabels['M0'] = 'm0'
xlabels['FWHM'] = 'fwhm ["]'

directory = '/mnt/scratch-lustre/njones/SURP2015/dflydata/'
subdir = '/darksub_flatfield/'
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
totalstats = []
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
                    shapearr = zeros((607,437))
                    time = Time(header['DATE'])
                    alt,az,x,y = getAltAz(shapearr,header,time,telescope)
                    midind = where((x == shapearr.shape[1]/2) & 
                                   (y == shapearr.shape[0]/2))
                    midalt = alt[midind]
                    X = airmass(midalt)
                    colour[serialno] = header['FILTNAM']
                    try:
                        datedict[serialno+'_'+expn].append(float(header[chosenkey]))
                    except (KeyError,NameError):
                        datedict[serialno+'_'+expn] = [float(header[chosenkey])]
                    try:
                        camdict[serialno].append(float(header[chosenkey]))
                    except (KeyError,NameError):
                        camdict[serialno] = [float(header[chosenkey])]
                    row = [key,str(time.jd),serialno,header['FILTNAM'],expn,
                           str(midalt[0]), str(X[0]),str(header['M0']),
                           str(header['FWHM']),str(date),str(header['SLOPE'])]
                    totalstats.append(row)
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

T = array(totalstats)

savetxt('stats/stats.txt',T,fmt='%s')

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

objects = T[:,0]
times = T[:,1].astype(float)
serials = T[:,2]
colours = T[:,3]
expns = T[:,4]
altitudes = T[:,5].astype(float)
airmasses = T[:,6].astype(float)
m0s = T[:,7].astype(float)
seeings = T[:,8].astype(float)

spi = where(objects == 'spi1_1')
dra = where(objects == 'PGM_1_2')

keys = {}
keys['spi'] = spi
keys['dra'] = dra
names = {}
names['spi'] = 'Spider'
names['dra'] = 'Draco'
for key in names.keys():

#SPIDER
    plt.figure()
    plt.plot(times[keys[key]],m0s[keys[key]],'.')
    plt.xlabel('Time [JD]')
    plt.ylabel('$m_0$')
    plt.title(names[key])
    plt.figure()
    plt.plot(times[keys[key]],seeings[keys[key]],'.')
    plt.xlabel('Time [JD]')
    plt.ylabel('FWHM')
    plt.title(names[key])
    plt.figure()
    plt.plot(m0s[keys[key]],seeings[keys[key]],'.')
    plt.xlabel('$m_0$')
    plt.ylabel('FWHM')
    plt.title(names[key])
