from numpy import *
from astropy import wcs
from astropy.io import fits
from cartesian import cartesian

def gridtowcs(fname,ind=0):
    f = fits.open(fname)
    header = f[ind].header
    data = f[ind].data
    f.close()
    w = wcs.WCS(header)
    x = arange(data.shape[0])
    y = arange(data.shape[1])
    coordpairs = cartesian([x,y])
    world = w.wcs_pix2world(coordpairs,0)
    return coordpairs, world

def pixtowcs(wcscoords,fname,ind):
    f = fits.open(fname)
    header = f[ind].header
    f.close()
    w = wcs.WCS(header)
    coords = w.wcs_pix2world(wcscoords,0)
    return coords

def wcstogrid(wcscoords,fname,ind):
    f = fits.open(fname)
    header = f[ind].header
    f.close()
    w = wcs.WCS(header)
    pixels = w.wcs_world2pix(wcscoords,0)
    return pixels
