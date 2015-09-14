"""
regrid - contains functions to do regriding and reshaping of image data

Requires the following modules: astropy, numpy, docopt, scipy, time,
                                collections
Requires the following files:   cartesian.py

Contains the following functions: regrid_lowmemory, regrid, nanlimits,
                                  valfilter, reshapeparams, reshape
"""

from astropy.io import fits
from astropy import wcs
from numpy import *
import docopt
import scipy.interpolate
import time
from cartesian import cartesian
from collections import Counter, OrderedDict

def regrid_lowmemory(sourceimage,targetimage,fillval = NAN,theader = 0):
    """
    A low memory version of regrid, this takes sourceimage and puts it onto
        targetimage grid, using wcs solutions for both. Grid points with 
        no info from sourceimage are set to fillval.

    sourceimage:    path to image to regrid - should already be convolved 
                    to appropriate resolution using resconvolve
    targetimage:    path to image whose grid is to be used in regridding
    fillval:        value to give to empty grid positions
                    (kwarg, default = NAN)
    theader:        specify a header containing WCS solution to use
                    (kwarg, default = 0)

    Returns array with targetimage dimensions.       

    """
    # Start timer
    start = time.time()
    # Load in source data and header information
    sdata,sheader = fits.getdata(sourceimage,header = True)
    # Create WCS object for the source image
    sourcewcs = wcs.WCS(sheader)
    # Create array of pixel indices in source image
    x = arange(sdata.shape[1])
    y = arange(sdata.shape[0])
    # Interpolate the image data over pixel indices
    interp = scipy.interpolate.RectBivariateSpline(y,x,sdata)
    # Load in target grid data
    if theader == 0:
        tdata,theader = fits.getdata(targetimage,header=True)
    if theader != 0:
        tdata = fits.getdata(targetimage)
    # Create WCS object for target grid
    targetwcs = wcs.WCS(theader)
    # Create all possible pairs of pixel coordinates in target grid
    coords = cartesian([arange(tdata.shape[1]),arange(tdata.shape[0])])
    # Extract x and y columns of pixel pairs
    xpixs = coords[:,0]
    ypixs= coords[:,1]
    # Convert target grid pixels to ra/dec 
    world = targetwcs.wcs_pix2world(coords,0)
    # Convert target grid ra/dec to source pixel coordinates
    dpix = sourcewcs.wcs_world2pix(world,0)
    # Extract x and y columns of converted pixel pairs
    xdpixs = dpix[:,0]
    ydpixs = dpix[:,1]
    # Find where target grid corresponds to actual source image data
    good = where((xdpixs >= min(x)) & (xdpixs <= max(x)) & (ydpixs >= min(y)) & (ydpixs <= max(y)))
    # Pick out only pixels with relevant image data
    xpixs = xpixs[good]
    ypixs = ypixs[good]
    xdpixs = xdpixs[good]
    ydpixs = ydpixs[good]
    # Create grid to fill up with source image regrid
    tofill = copy(tdata)
    tofill[:] = fillval
    # Loop over relevant pixels and fill up source image regrid
    for i in range(len(ypixs)):
        ypix = ypixs[i]
        xpix = xpixs[i]
        xdpix = xdpixs[i]
        ydpix = ydpixs[i]
        tofill[ypix,xpix] = interp(ydpix,xdpix)
    # End timer
    end = time.time()
    # Print time to run
    print 'Regridded in ',(end-start)/60.,' min'
    return tofill,tdata


def regrid(sourceimage,targetimage,fillval = NAN,theader = 0,tpix = []):
    """
    This takes sourceimage and puts it onto targetimage grid, using wcs 
        solutions for both. Grid points with no info from sourceimage are 
        set to fillval. Requires scipy0.14.0

    sourceimage:    path to image to regrid - should already be convolved 
                    to appropriate resolution using resconvolve
    targetimage:    path to image whose grid is to be used in regridding
    fillval:        value to give to empty grid positions
                    (kwarg, default = NAN)
    theader:        specify a header containing WCS solution to use
                    (kwarg, default = 0)

    Returns array with targetimage dimensions.      
    """
    # Start timer
    start = time.time()
    # Load in source data and header information
    sdata,sheader = fits.getdata(sourceimage,header = True)
    # Create WCS object for the source image
    sourcewcs = wcs.WCS(sheader)
    # Create array of pixel indices in source image
    x = arange(sdata.shape[1])
    y = arange(sdata.shape[0])
    # Interpolate the image data over pixel indices
    interp = scipy.interpolate.RectBivariateSpline(y,x,sdata)
    # Load in target grid data
    if theader == 0 or tpix == []:
        tdata,theader = fits.getdata(targetimage,header=True)
        # Create all possible pairs of pixel coordinates in target grid
        coords = cartesian([arange(tdata.shape[1]),arange(tdata.shape[0])])
    elif theader != 0 and tpix != []:
        assert len(tpix) == 4
        dx,ux,dy,uy = tpix
        # Create all possible pairs of pixel coordinates in target grid
        tdata = zeros((uy-dy,ux-dx))
        coords = cartesian([arange(dx,ux),arange(dy,uy)])
    # Extract x and y columns of pixel pairs
    xpixs = coords[:,0]
    ypixs= coords[:,1]
    # Create WCS object for target grid
    targetwcs = wcs.WCS(theader)
    # Convert target grid pixels to ra/dec 
    world = targetwcs.wcs_pix2world(coords,0)
    # Convert target grid ra/dec to source pixel coordinates
    dpix = sourcewcs.wcs_world2pix(world,0)
    # Extract x and y columns of converted pixel pairs
    xdpixs = dpix[:,0]
    ydpixs = dpix[:,1]
    # Find where target grid corresponds to actual source image data
    good = where((xdpixs >= min(x)) & (xdpixs <= max(x)) & 
                 (ydpixs >= min(y)) & (ydpixs <= max(y)))
    # Pick out only pixels with relevant image data
    xpixs = xpixs[good]
    ypixs = ypixs[good]
    xdpixs = xdpixs[good]
    ydpixs = ydpixs[good]
    # Create grid to fill up with source image regrid
    tofill = copy(tdata)
    tofill[:] = fillval
    # Choose indices of array positions to be changed
    inds = (ypixs,xpixs)
    tofill[inds] = interp(ydpixs,xdpixs,grid=False) # needs scipy 0.14.0
    # End timer
    end = time.time()
    # Print time to run
    print 'Regridded in ',(end-start)/60.,' min'
    return tofill,tdata

def vallimits(arr,limval = NAN):
    """
    For 1D arr, determine the index location of the first
        and last non-NAN value.

    arr:    1D input array

    Returns 2 indices

    """
    # Find where arr is non-NAN
    if isnan(limval) == True:
        vallocs = where(isnan(arr) == False)
    else:
        vallocs = where(arr != limval)
    # Extract the indices corresponding to the maximum and
    # minimum non-NANs. If none exist, set them to NAN
    try:
        uplim = max(vallocs[0])
    except ValueError:
        uplim = NAN
    try:
        downlim = min(vallocs[0])
    except ValueError:
        downlim = NAN
    return downlim,uplim

def valfilter(ls,minmax):
    """
    Filters a list by finding an extremal value that occurs more than once

    ls:     list to filter
    minmax: specify whether minimum or maximum value required

    Returns a float

    """
    # Find how often each values occur in ls
    count = Counter(ls)
    # Remove keys that occur only once
    keys = count.keys()
    for key in keys:
        if count[key] == 1:
            del count[key]
    keys = count.keys()
    # Return min or max as specified
    if minmax == 'min':
        return min(keys)
    if minmax == 'max':
        return max(keys)

def reshapeparams(data,limval = NAN):
    """
    Finds slicing limit of 2D array data

    data:   2D array to be sliced

    Returns slicing limits, 4 indices
    """
    # Prepare list of potential indices
    rowups = []
    colups = []
    rowdowns = []
    coldowns = []
    # Search rows for limits
    for i in range(len(data)):
        rowdown,rowup = vallimits(data[i],limval = limval)
        if not isnan(rowup):
            rowups.append(rowup)
        if not isnan(rowdown):
            rowdowns.append(rowdown)
    # Search columns for limits
    for i in range(len(data[0])):
        coldown,colup = vallimits(data[:,i],limval = limval)
        if not isnan(colup):
            colups.append(colup)
        if not isnan(coldown):
            coldowns.append(coldown)
    # Of list of potential indices, find min and max
    rowup = valfilter(rowups,'min')
    rowdo = valfilter(rowdowns,'max')
    colup = valfilter(colups,'min')
    coldo = valfilter(coldowns,'max')
    return rowup,rowdo,colup,coldo

def reshape(data,shapebyarr,ress,limval = NAN):
    """
    Reshapes 2D data according to 2D shapebyarr, then slices ress rows
        and columns off of each side.

    data:       2D array to be reshaped
    shapebyarr: 2D array to serve as instructions for reshaping
    ress:       number of pixels to remove from each side

    Returns a 2D array
    """
    # Find the slicng limits for the array
    rowup,rowdo,colup,coldo = reshapeparams(shapebyarr,limval = limval)
    # Slice array by NAN limits and additional limits
    colslice = data[coldo+ress:colup-ress]
    rowslice = colslice.T
    rowslice = rowslice[rowdo+ress:rowup-ress]
    newdata = rowslice.T
    return newdata
