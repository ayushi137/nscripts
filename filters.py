from astropy.io import fits
import os

directory = '/mnt/scratch-lustre/njones/SURP2015/dflydata/'
subdir = '/darksub_flatfield/'

objectdates = {}
objectdates['spi1_1'] = ['2014-11-19','2014-11-22','2014-11-23','2014-11-24']

cams = {}
cams['83F010687'] = 'SloanG'
cams['83F010730'] = 'SloanG'
cams['83F010783'] = 'SloanG'
cams['83F010820'] = 'SloanG'
cams['83F010827'] = 'SloanG'

cams['83F011129'] = 'SloanR'
cams['83F010692'] = 'SloanR'
cams['83F010784'] = 'SloanR'
cams['83F010826'] = 'SloanR'

for objects in objectdates.keys():
    dates = objectdates[objects]
    for date in dates:
        di = directory+objects+subdir+date+'/'
        files = os.listdir(di)
        for f in files:
            if os.path.isfile(di+f):
                print di+f
                data,header = fits.getdata(di+f,header=True)
                serialno = f[0:9]
                header['FILTNAM'] = cams[serialno]
                fits.writeto(di+f,data,header,clobber = True)

