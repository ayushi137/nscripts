"""
maskdata - contains a function to mask image data

Requires the following modules: os, numpy, docopt, astropy
                                matplotlib

Contains the following functions: maskdata

Usage:
maskdata [-h]

Options:
    -h, --help         Show this screen

"""

########################### IMPORT BASE PACKAGES ############################

import os
from numpy import *
import docopt
from astropy.io import fits
import matplotlib.pyplot as plt

docopt.docopt(__doc__)

########################### FUNCTIONS ############################

def maskdata(data,starmap,cutoff,masksave = 0,outfile=0,header=0):
    """
    Masks data according to starmap

    data:           image to be masked
    starmap:        image to use as mask
    cutoff:         value to use as cutoff in starmap, below which pixels are unmasked
    masksave:       name of file to save mask into if value is zero,
                    do not save file
                    (kwarg, default = 0)
    outfile:        name of outfile, if masked image is to be 
                    saved - if value is zero, do not save file
                    (kwarg, default = 0)
    header:         header information, if masked image is to
                    be saved - if value is zero, do not save file
                    (kwarg, default = 0)

    Returns masked data.
    """
    # Find indices above and below cutoff
    above = where(starmap >= cutoff)
    below = where(starmap < cutoff)
    # Create mask from starmap
    mask = copy(starmap)
    mask[above] = 0
    mask[below] = 1
    # Save mask, if required
    if masksave != 0 and header != 0:
        fits.writeto(masksave,mask,header,clobber=True)
    # Mask the image
    maskeddata = mask*data
    # Save the file, if necessary
    if outfile != 0 and header != 0:
        header['MASKCUT'] = cutoff
        fits.writeto(outfile,maskeddata,header,clobber=True)
        return maskeddata,header
    elif outfile == 0 or header == 0:
        return maskeddata,0
