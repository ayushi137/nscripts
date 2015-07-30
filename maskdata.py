import os
from numpy import *
from astropy.io import fits
import matplotlib.pyplot as plt

def maskdata(data,starmap,cutoff,outfile=0,header=0):
    above = where(starmap >= cutoff)
    below = where(starmap < cutoff)
    mask = copy(starmap)
    mask[above] = 0
    mask[below] = 1
    maskeddata = mask*data
    if outfile != 0 and header != 0:
        header['MASKCUT'] = cutoff
        fits.writeto(outfile,maskeddata,header,clobber=True)
        return maskeddata,header
    elif outfile == 0 or header == 0:
        return maskeddata
