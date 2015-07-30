from astropy.io import fits
from astropy import wcs
from numpy import *
import scipy.interpolate
import time
from cartesian import cartesian
from collections import Counter, OrderedDict

def regridoldschool(sourceimage,targetimage,fillval = NAN):
    """
    sourceimage = path to image to regrid - should already be convolved to appropriate 
                  resolution using resconvolve
    targetimage = path to image whose grid is to be used in regridding
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
    tdata,theader = fits.getdata(targetimage,header=True)
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

def nanlimits(array):
    vallocs = where(isnan(array) == False)
    try:
        uplim = max(vallocs[0])
    except ValueError:
        uplim = NAN
    try:
        downlim = min(vallocs[0])
    except ValueError:
        downlim = NAN
    return downlim,uplim

def regrid(sourceimage,targetimage,fillval = NAN):
    """
    ****************REQUIRES SCIPY0.14.0 AT MIN **************************
    sourceimage = path to image to regrid - should already be convolved to 
                  appropriate resolution using resconvolve
    targetimage = path to image whose grid is to be used in regridding
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
    tdata,theader = fits.getdata(targetimage,header=True)
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

def valfilter(ls,minmax):
    count = Counter(ls)
    ocount = OrderedDict(sorted(count.items()))
    keys = ocount.keys()
    values = ocount.values()
    for key in keys:
        if ocount[key] == 1:
            del ocount[key]
    keys = ocount.keys()
    if minmax == 'min':
        return min(keys)
    if minmax == 'max':
        return max(keys)

def reshapeparams(data):
    # search rows
    rowups = []
    colups = []
    rowdowns = []
    coldowns = []
    for i in range(len(data)):
        rowdown,rowup = nanlimits(data[i])
        if not isnan(rowup):
            rowups.append(rowup)
        if not isnan(rowdown):
            rowdowns.append(rowdown)
    for i in range(len(data[0])):
        coldown,colup = nanlimits(data[:,i])
        if not isnan(colup):
            colups.append(colup)
        if not isnan(coldown):
            coldowns.append(coldown)
    rowup = valfilter(rowups,'min')
    rowdo = valfilter(rowdowns,'max')
    colup = valfilter(colups,'min')
    coldo = valfilter(coldowns,'max')
    return rowup,rowdo,colup,coldo

def reshape(data,shapebyarr,ress):
    rowup,rowdo,colup,coldo = reshapeparams(shapebyarr)
    colslice = data[coldo+ress:colup-ress]
    rowslice = colslice.T
    rowslice = rowslice[rowdo+ress:rowup-ress]
    newdata = rowslice.T
    return newdata
