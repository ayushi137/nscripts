from regrid import regrid
from astropy.io import fits
import matplotlib.pyplot as plt
from numpy import *
plt.ion()

def coadd(flist,standard,outname=0):
    sdata,sheader = fits.getdata(standard,header = True)
    cube = zeros((len(flist),sdata.shape[0],sdata.shape[1]))
    cube[0] = sdata             
    for i in range(len(flist)):
        print i+1, ' of ', len(flist)
        re,target = regrid(flist[i],standard,fillval=0)
        plt.figure()
        plt.imshow(re)
        cube[i] = re
    cadd = cube.mean(axis = 0)
    plt.figure()
    plt.imshow(cadd)
    sheader['COADD'] = 'TRUE'
    if outname != 0:
        fits.writeto(outname,cadd,sheader,clobber=True)
    elif outname == 0:
        return cadd
