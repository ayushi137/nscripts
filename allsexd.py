import os
from numpy import *
from callastrometry import callastrometry, scrubwcsheader

def sexcalls(f,di,odi,cdi,bdi):
    fname = f.split('.fits')[0]
    objsexcall = 'sex -PARAMETERS_NAME photo.param -CATALOG_NAME '+cdi+fname+'.cat'+' -CHECKIMAGE_TYPE OBJECTS, BACKGROUND -CHECKIMAGE_NAME '+odi+fname+'_objects.fits '+di+f+', '+bdi+fname+'_background.fits '+di+f
    os.system(objsexcall)

if __name__ == '__main__':
    directory = '/mnt/gfsproject/naiad/njones/moddragonflydata/'
    subdir = '/NCalibrated/'
    objdir = '/objects/'
    bakdir = '/background/'
    catdir = '/catalogue/'
    dirlist = [objdir,bakdir,catdir]
    for d in dirlist:
        if os.path.isdir(directory+'PGM_1_2'+d) == False:
            os.system('mkdir '+directory+'PGM_1_2'+d)
        if os.path.isdir(directory+'spi1_1'+d) == False:
            os.system('mkdir '+directory+'spi1_1'+d)
    
    objectdates = {}
    objectdates['PGM_1_2'] = ['2013-09-23','2013-09-24','2013-09-25',
                              '2013-09-26','2013-09-27']
    objectdates['spi1_1'] = ['2014-10-31','2014-11-01','2014-11-19',
                             '2014-11-22','2014-11-23','2014-11-24']
    
    keys = objectdates.keys()
    #keys = ['PGM_1_2']
    for key in keys:
        dates = objectdates[key]
        #dates = ['2013-09-23']
        for date in dates:
            di = directory+key+subdir+date+'/'
            files = os.listdir(di)
            odi = directory+key+objdir+date+'/'
            bdi = directory+key+bakdir+date+'/'
            cdi = directory+key+catdir+date+'/'
            dis = [odi,bdi,cdi]
            for d in dis:
                if os.path.isdir(d) == False:
                    os.system('mkdir '+d)
            files.sort()
            #files = ['83F010826_17_light_ds_ff.fits']
            for f in files:
                if '.wcs' not in f:
                    print f
                    sexcalls(f,di,odi,cdi,bdi)
                
