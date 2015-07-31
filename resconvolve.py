"""
resconvolve - contains a function to use convolution to change the 
    resolution of an image

Requires the following modules: docopt, scipy, numpy, os, astropy

Contains the following functions: resconvolve

Usage:
resconvolve [-h]

Options:
    -h, --help      Show this screen

"""

########################## IMPORT PACKAGES ###########################

import docopt
import scipy.signal
import scipy.interpolate
from numpy import *
import os
from astropy.io import fits

docopt.docopt(__doc__)

########################## FUNCTIONS ###########################

def resconvolve(data,currentres,desiredres,pixscale = 2.85,
                outfile=0,header=0):
    """
    Convolves data from currentres to desiredres.

    data:           image to change resolution of
    currentres:     current resolution of the image in arcseconds
    desiredres:     desired resolution to change to in arcseconds
    pixscale:       either float or 2-element list, the pixel scale
                    in arcseconds/pixel - if a list, first element    
                    is x-dimension pixel scale and second is y
                    (kwarg, default = 2.85)
    outfile:        name of outfile, if convolved image is to be 
                    saved - if value is zero, do not save file
                    (kwarg, default = 0)
    header:         header information, if convolved image is to
                    be saved - if value is zero, do not save file
                    (kwarg, default = 0)

    Returns convolved data

    """
    # Factor to convert sigma used in scipy.signal.gaussian to FWHM
    s2f = 2*sqrt(2*log(2))
    # Convert resolutions to sigma for scipy.signal.gaussian
    currentres /= s2f
    desiredres /= s2f
    # Check that convolution to change resolution is sensible   
    if  desiredres < currentres:
        print 'Convolution failure'
        print 'current resolution: ',currentres
        print 'desired resolution: ',desiredres
        return data,0
    elif desiredres > currentres:
        # Find x-dimension
        gridx = data.shape[0]
        gridy = data.shape[1]
        # Determine the size of the kernel in arcseconds
        kernelres = sqrt(desiredres**2-currentres**2)
        # Convert kernel to pixels and create normalized 2D Gaussian
        if isinstance(pixscale,float):
            kernelsize = kernelres/pixscale
            pre = 1./(kernelsize*sqrt(2*pi))
            kernelx = pre*scipy.signal.gaussian(gridx,kernelsize)
            kernely = pre*scipy.signal.gaussian(gridy,kernelsize)
        elif isinstance(pixscale,(list,ndarray)):
            kernelxsize = kernelres/pixscale[0]
            kernelysize = kernelres/pixscale[1]
            prex = 1./(kernelxsize*sqrt(2*pi))
            prey = 1./(kernelysize*sqrt(2*pi))
            kernelx = prex*scipy.signal.gaussian(gridx,kernelxsize)
            kernely = prey*scipy.signal.gaussian(gridy,kernelysize)            
        kernel = outer(kernelx,kernely)
        # Convolve image data and kernel using fft
        convolved = scipy.signal.fftconvolve(data,kernel,mode = 'same')
        # Update header and save file if necessary, then return data
        if outfile != 0 and header != 0:
            if isinstance(pixscale,float):
                header['CONVKER'] = (kernelsize,'sigma of convolution gaussian in dfly pix')
            elif isinstance(pixscale,(list,ndarray)):
                header['CONVKERX'] = (kernelxsize,'sigma of x-convolution gaussian in dfly pix')
                header['CONVKERY'] = (kernelysize,'sigma of y-convolution gaussian in dfly pix')
            fits.writeto(outfile,convolved,header,clobber = True)
            return convolved,header
        elif header == 0 or outfile == 0:
            return convolved,header
