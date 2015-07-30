import os
from astropy import wcs
from astropy.io import fits

fileextensions = ['.axy','.corr','-indx.xyls','.match','.new','-objs.png',
                  '.rdls','.solved']
def callastrometry(fname,generate=False):
    keys = ['CTYPE1','CTYPE2','CRPIX1','CRVAL1','CRPIX2','CRVAL2','CD1_1',
            'CD1_2','CD2_1','CD2_2','EQUINOX']
    header = fits.getheader(fname)
    if set(keys) < set(header.keys()) and generate == False:
        print 'Astrometry already done'
    else:
        command = 'solve-field --no-fits2fits --use-sextractor --cpulimit 20 {0}'.format(fname)
        os.system(command)
        trimfname = fname.split('.fits')[0]
        print 'TRIMFNAME = '+trimfname
        for ext in fileextensions:
            os.system('rm -f '+trimfname+ext)
        lightwcs = fits.getheader(trimfname+'.wcs')
        w = wcs.WCS(lightwcs)
        scales = wcs.utils.proj_plane_pixel_scales(w)*3600
        data,light = fits.getdata(fname,header=True)
        light['PSCALX']=scales[0]
        light['PSCALY']=scales[1]
        light['CTYPE1']=lightwcs['CTYPE1']
        light['CTYPE2']=lightwcs['CTYPE2']
        light['CRPIX1']=lightwcs['CRPIX1']
        light['CRVAL1']=lightwcs['CRVAL1']
        light['CRPIX2']=lightwcs['CRPIX2']
        light['CRVAL2']=lightwcs['CRVAL2']
        light['CD1_1']=lightwcs['CD1_1']
        light['CD1_2']=lightwcs['CD1_2']
        light['CD2_1']=lightwcs['CD2_1']
        light['CD2_2']=lightwcs['CD2_2']
        light['EQUINOX']=lightwcs['EQUINOX']
        fits.writeto(fname,data,light,clobber=True)
        os.system('rm -f '+trimfname+'.wcs')

def scrubwcsheader(fname,ind=0):
    keys = ['CTYPE1','CTYPE2','CRPIX1','CRVAL1','CRPIX2','CRVAL2','CD1_1',
            'CD1_2','CD2_1','CD2_2','EQUINOX','PSCALX','PSCALY']
    light = fits.open(fname)
    header = light[ind].header
    for key in keys:
        try:
            del header[key]
        except (ValueError,KeyError):
            continue
    data = light[ind].data
    fits.writeto(fname,data,header,clobber = True)
