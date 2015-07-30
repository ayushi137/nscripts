import scipy.signal
import scipy.interpolate
from numpy import *
import os
from astropy.io import fits

def resconvolve(data,currentres,desiredres,pixscale = 2.85,outfile=0,header=0):
    if  desiredres < currentres:
        print 'Convolution failure'
        print 'current resolution: ',currentres
        print 'desired resolution: ',desiredres
        return data,0
    elif desiredres > currentres:
        kernelres = sqrt(desiredres**2-currentres**2)
        kernelsize = kernelres/pixscale
        gridx = data.shape[0]
        gridy = data.shape[1]
        pre = 1./(kernelsize*sqrt(2*pi))
        kernelx = pre*scipy.signal.gaussian(gridx,kernelsize)
        kernely = pre*scipy.signal.gaussian(gridy,kernelsize)
        kernel = outer(kernelx,kernely)
        convolved = scipy.signal.fftconvolve(data,kernel,mode = 'same')
        if outfile != 0 and header != 0:
            header['CONVKER'] = (kernelsize,'sigma of convolution gaussian in dragonfly pixels')
            fits.writeto(outfile,convolved,header,clobber = True)
            return convolved,header
        elif header == 0 or outfile == 0:
            return convolved
